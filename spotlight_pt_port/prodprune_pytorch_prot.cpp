#define STATICDATA 1
#define NOJSON 1
#define NOINTERFACE 1
#define USE_TORCH 1
//#define FFCHARGES 1 //Uncomment to disable machine learned charges.
//#define SERIAL 1 //Enable to run serial job

//#ifndef SERIAL
//#include <omp.h>
//#endif

#include "graph/Molecule.hpp"

static std::string PTPATH="../";

//ML part (start)
static bool RANDOM_ATOM=false;
#ifndef SERIAL
  static const int MAX_THREADS_OMP=omp_get_max_threads(); //Change here if needed
#endif

/* A simple 2-layer MLP model */
class SimpleModelImpl : public mytorch::TorchModuleWrapper
{
public:
    SimpleModelImpl(int64_t D_in, int64_t H) : mytorch::TorchModuleWrapper()
    {
        // Construct and register two Linear submodules.
        fc1 = register_module("fc1", torch::nn::Linear(D_in, H));
        fc2 = register_module("fc2", torch::nn::Linear(H, 1));
    }

    // Implement the Net's algorithm.
    torch::Tensor forward(mytorch::TorchTensorList xv) override
    {
        // Use one of many tensor manipulation functions.
        torch::Tensor x = torch::relu(fc1->forward(xv[0]));
        x = fc2->forward(x);
        return x;
    }

    // Use one of many "standard library" modules.
    torch::nn::Linear fc1{nullptr}, fc2{nullptr};
};
TORCH_MODULE(SimpleModel);

class RNNModelImpl : public mytorch::TorchModuleWrapper
{
  torch::nn::GRU rnn{nullptr};
  torch::nn::Linear linear{nullptr};
  torch::Tensor memory;
  int H;
public:
  RNNModelImpl(int64_t D_in, int64_t H) : mytorch::TorchModuleWrapper()
  {
    this->H=H;
    rnn = register_module("rnn1",torch::nn::GRU(torch::nn::GRUOptions(D_in,H).num_layers(1).bidirectional(false).batch_first(true)));
    linear=register_module("linear1",torch::nn::Linear(H,1));
    this->resetModel();
  }

  torch::Tensor forward(mytorch::TorchTensorList xv) override
  {
    //std::cout << "Hidden Size: "<<mytorch::getShape(memory)<<"\n";
    torch::Tensor xup=xv[0].unsqueeze(0);
    auto rnnout=rnn->forward(xup,memory);
    memory=std::get<1>(rnnout);
    //memory=rnnout.state;
    torch::Tensor output= linear->forward(std::get<0>(rnnout)).squeeze(0);
    //std::cout <<"Output Size: "<<mytorch::getShape(output)<<"\n";
    return output;
    //return linear->forward(rnnout.output);
  }

  inline void resetModel() {memory=torch::zeros({1,1,H});}
};
TORCH_MODULE(RNNModel);

static const int BLOCK_SIZE=50;
static std::vector<torch::Tensor> builtmols;

static ForceField ff; 
static float randomchance=0.0;

static Atom* internalGetPolicyAtom(const ForceField& ff,const std::vector<std::string>& atlst,Atom* src,const std::vector<Atom*>& nhd,int devid=0)
{
  torch::NoGradGuard no_grad;
  std::vector<Atom*> chs=ff.toAtomList(atlst);
  torch::Tensor allowed_atomfeats=featurizeAtomTypes(*atom_featurizer,chs);
  #ifdef SERIAL
  torch::Tensor sing=molenvmt->getMoleculeFeatures()[molenvmt->getMolecule()->indexOf(src)].repeat({allowed_atomfeats.size(0),1}); //atom_featurizer->featurize(src,nhd).repeat({allowed_atomfeats.size(0),1});
  #else
  torch::Tensor sing=molenvmt[devid]->getMoleculeFeatures()[molenvmt[devid]->getMolecule()->indexOf(src)].repeat({allowed_atomfeats.size(0),1}); //atom_featurizer->featurize(src,nhd).repeat({allowed_atomfeats.size(0),1});
  #endif
  allowed_atomfeats = torch::cat({allowed_atomfeats,sing},1);
  int retidx = 0;
  if(RANDOM_ATOM && randomchance>0 && RNG::throwarandompoint(0,1)<randomchance) retidx=(int)(RNG::throwarandompoint(0,1)*chs.size());
  #ifdef SERIAL
  else retidx=mypolicy->select_action(allowed_atomfeats);
  mypolicy->save_reward(0.0f);
  #else
  else retidx=mypolicy[devid]->select_action(allowed_atomfeats);
  mypolicy[devid]->save_reward(0.0f);
  #endif
  return chs[retidx];
}

static int MOL_REG_CONST=125,MOL_REG_THRESH=5;
static float VARIANCE_REG_CONST=0.01f,REG_SCALE=12.5f,VARIANCE_REG_SCALE=1.250f;
static torch::Tensor getAtomTypeWeights(torch::Tensor t)
{
	torch::Tensor cnts=torch::sum(t,0);
	float norm=torch::sum(t).item<float>();
	return cnts/norm;
}

//ML part (end)

struct InputParameters
{
  std::string fnam="",format="pdb",progname="dnv_protgen",flagcodes="0000";
  int num=25,osc=5,seednum=30,step=1;
  std::vector<int> size={15};
  std::string source_folder="data";
  std::string ff_file=source_folder+"/final_ff_parameters.ffin",cats_file=source_folder+"/categories.data",charges_file=source_folder+"/charges.data",defs_file=source_folder+"/definitions.data";
  std::string bondsfile=source_folder+"/itps/bondtypes.itp",anglesfile=source_folder+"/itps/angletypes.itp",dihedralsfile=source_folder+"/itps/dihedraltypes.itp",impropersfile=source_folder+"/itps/impropertypes.itp";
  std::string protein_ff=source_folder+"/protff.ffin",resid_ff=source_folder+"/protein.rtp";
  std::string extra_atypes="",extra_rules="";
  std::string extra_bonds="",extra_angles="",extra_dihs="";
  double dielectric=30,restrain_distance=0.8,spread=1.2;
  float prunefrac=0.67;
  std::string prunemolfile="prunesrc.pdb";
  bool prune=false;

  inline void updateFFData()
  {
    ff_file=source_folder+"/final_ff_parameters.ffin";
    cats_file=source_folder+"/categories.data";
    charges_file=source_folder+"/charges.data";
    defs_file=source_folder+"/definitions.data";
    bondsfile=source_folder+"/itps/bondtypes.itp";
    anglesfile=source_folder+"/itps/angletypes.itp";
    dihedralsfile=source_folder+"/itps/dihedraltypes.itp";
    impropersfile=source_folder+"/itps/impropertypes.itp";
    protein_ff=source_folder+"/protff.ffin";
    resid_ff=source_folder+"/protein.rtp";
  }
};
static std::vector<int> getResidueNumbersFromIdentifiers(const std::vector<std::string>& residIDs,int step=1)
{
  std::vector<int> ret;
  for(std::string s : residIDs)
  {
    if(stringfx::indexOf('-',s)!=std::string::npos)
    {
      auto parts=stringfx::split(s,'-');
      int start=std::stoi(parts[0]),end=std::stoi(parts[1]);
      for(int ie=start;ie<=end;ie+=step) ret.push_back(ie);
    }
    else ret.push_back(std::stoi(s));
  }
  return ret;
}
int main(int argc,char** argv)
{
  torch::NoGradGuard no_grad;
  std::string filename=argv[1];
  std::ifstream loadf; loadf.open(filename);
  InputParameters INP;
  std::string line;
  std::vector<std::string> residIDs;
  std::vector<std::string> seednames;

  //All options defaults
  char s1='0',s2='0',s3='0',s4='0';
  USE_PHYSICS=false;
  QUIET=true;
  SYNTHGET=false;
  IPADDR="127.0.0.1";
  PORT=5555;
  PRUNEKEEPHS=false;
  #ifdef USE_TORCH
  WEIGHTED_ATOM_SELECTOR_FUNCTION=internalGetPolicyAtom;
  #endif
  //SYNTHPENALTY=350; //Strong - Don't accept without //Default - Set very high to get exclusively synthesizable molecules (mostly)

  std::string sizes="15";


  while(getline(loadf,line))
  {
    std::vector<std::string> tokens=stringfx::split(line,'=');
    std::string idtoken=stringfx::trim(tokens[0]);
    if(idtoken=="Protein") INP.fnam=stringfx::trim(tokens[1]);
    else if(idtoken=="ProteinFormat") INP.format=stringfx::trim(tokens[1]);
    else if(idtoken=="Dielectric") INP.dielectric=std::stod(stringfx::trim(tokens[1]));
    else if(idtoken=="RestrainDistance") INP.restrain_distance=std::stod(stringfx::trim(tokens[1]));
    else if(idtoken=="SeedCount") INP.seednum=std::stoi(stringfx::trim(tokens[1]));
    else if(idtoken=="Step") INP.step=std::stoi(stringfx::trim(tokens[1]));
    else if(idtoken=="Sizes" || idtoken=="Sizes") sizes=stringfx::trim(tokens[1]);
    else if(idtoken=="MolCount") INP.num=std::stoi(stringfx::trim(tokens[1]));
    else if(idtoken=="SourceFolder") INP.source_folder=stringfx::trim(tokens[1]);
    else if(idtoken=="Oscillations") INP.osc=std::stoi(stringfx::trim(tokens[1]));
    else if(idtoken=="ProgramName") INP.progname=stringfx::trim(tokens[1]);
    else if(idtoken=="Strategy_Restrain") s1=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
    else if(idtoken=="Strategy_Superseed") s2=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
    else if(idtoken=="Strategy_Reseed") s3=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
    else if(idtoken=="Strategy_Deseed") s4=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
    else if(idtoken=="Residues") residIDs=stringfx::split(tokens[1],',');
    else if(idtoken=="PruneMol") {INP.prunemolfile=stringfx::trim(tokens[1]); INP.prune=true;}
    else if(idtoken=="PruneKeepFrac") INP.prunefrac=std::stof(stringfx::trim(tokens[1]));
    else if(idtoken=="SynthCheck") SYNTHGET=(stringfx::trim(tokens[1])=="Yes");
    else if(idtoken=="SynthStrength") SYNTHPENALTY=std::stoi(stringfx::trim(tokens[1]));
    else if(idtoken=="Optimize") USE_PHYSICS=(stringfx::trim(tokens[1])=="Yes");
    else if(idtoken=="Quiet") QUIET=(stringfx::trim(tokens[1])=="Yes");
    else if(idtoken=="Charges") FFCHARGES=(stringfx::trim(tokens[1])=="FF");
    else if(idtoken=="IP") IPADDR=stringfx::trim(tokens[1]);
    else if(idtoken=="PORT") PORT=std::stoi(tokens[1]);
    else if(idtoken=="SeedSpread" || idtoken=="Spread") INP.spread=std::stod(tokens[1]);
    else if(idtoken=="Seeds")
    {
      if(stringfx::trim(tokens[1])=="All") seednames={"CT1","CT2","CT3","CA","NX","CAS","NX","CPT","NY","CY","CT2x","CT1x","CP2","CE1A","CE2A","PL","SL","C3","OH1","OC"};
      else seednames=stringfx::split(tokens[1],',');
    }
    else if(idtoken=="ExtraAtoms") INP.extra_atypes+=","+stringfx::trim(tokens[1]);
    else if(idtoken=="ExtraBonds") INP.extra_bonds+=","+stringfx::trim(tokens[1]);
    else if(idtoken=="ExtraAngles") INP.extra_angles+=","+stringfx::trim(tokens[1]);
    else if(idtoken=="ExtraDihs" || idtoken=="ExtraDihedrals") INP.extra_dihs+=","+stringfx::trim(tokens[1]);
    else if(idtoken=="ExtraRules" || idtoken=="ExtraDefs" || idtoken=="ExtraDefinitions") INP.extra_rules+=","+stringfx::trim(tokens[1]);
  }
  INP.size=getResidueNumbersFromIdentifiers(stringfx::split(sizes,','),INP.step);
  std::string strategy(1,s1); strategy+=s2; strategy+=s3; strategy+=s4;
  INP.flagcodes=strategy;
  INP.updateFFData();
  int tempv;

 DIEL/=(INP.dielectric/3); //Dielectric constant is 30
 RESTRAINDISTANCE=INP.restrain_distance; //Restrain distance if enabled
 quickgeom::Cuboid* region=nullptr;

 //std::string fnam="",format="pdb",progname=(argc>7)?argv[7]:"dnv_protgen",flagcodes=(argc>8)?argv[8]:"0000";
 //int size=(argc>3)?std::stoi(argv[3]):15,num=(argc>4)?std::stoi(argv[4]):25,osc=(argc>5)?std::stoi(argv[5]):5,seednum=(argc>6)?std::stoi(argv[6]):30,tempv;
 //if(argc>1) fnam=argv[1];
 //if(argc>2) format=argv[2];
 //else {cout << "Enter filename: "; cin >>fnam;}
 cout << "Under: '"<<INP.progname<<"': Generating "<<INP.num<<" drugs of size="<<INP.size<<" using "<<INP.seednum<<" seeds, accumulating the best for every "<<INP.osc<<" oscillations.\n";
 cout << "Using the protein file: "<<INP.fnam<<"\tin the "<<INP.format<<" format.\n";
 if(INP.flagcodes!="0000") cout << "Running with: "<<((INP.flagcodes[0]!='0')?"restraint "+to_string(RESTRAINDISTANCE)+"A ":"")<<((INP.flagcodes[1]!='0')?"superseed ":"")<<((INP.flagcodes[2]!='0')?"reseed ":"")<<((INP.flagcodes[2]!='0')?"deseed ":"")<<"\n";
 else std::cout << "NO restraint, reseed, deseed, or superseed\n";
 std::cout << "Using dielectric: " << DIEL << "\n";
 std::vector<int> hsres=getResidueNumbersFromIdentifiers(residIDs);
 cout << "Using residue IDs: "<<hsres<<" with seeds: "<<seednames<<"\n";
 if(INP.prune)
 {
   cout << "PRUNE mode\n";
   cout << "Prune Keep Fraction: "<<INP.prunefrac<<"\n";
 }

 ForceField ff(INP.ff_file,INP.cats_file);
 ff.loadBondParameters(INP.bondsfile,INP.anglesfile,INP.dihedralsfile,INP.impropersfile);
 std::string extra_data="";
 if(INP.extra_atypes!="")
 {
   extra_data+="\tAtom types: "+INP.extra_atypes.substr(1)+"\n";
   for(const std::string fnm : stringfx::split(INP.extra_atypes.substr(1),',')) {std::cout <<"Loading: "<< fnm << "\n"; ff.loadAdditionalAtomTypes(stringfx::trim(fnm));}
 }
 if(INP.extra_bonds!="")
 {
   extra_data+="\tBond Data: "+INP.extra_bonds.substr(1)+"\n";
   for(const std::string fnm : stringfx::split(INP.extra_bonds.substr(1),',')) {std::cout <<"Loading: "<< fnm << "\n"; ff.loadLengths(stringfx::trim(fnm),false);}
 }
 if(INP.extra_angles!="")
 {
   extra_data+="\tAngle Data: "+INP.extra_angles.substr(1)+"\n";
   for(const std::string fnm : stringfx::split(INP.extra_angles.substr(1),',')) {std::cout <<"Loading: "<< fnm << "\n"; ff.loadAngles(stringfx::trim(fnm),false);}
 }
 if(INP.extra_dihs!="")
 {
   extra_data+="\tDihedral Data: "+INP.extra_dihs.substr(1)+"\n";
   for(const std::string fnm : stringfx::split(INP.extra_dihs.substr(1),',')) {std::cout <<"Loading: "<< fnm << "\n"; ff.loadDihedrals(stringfx::trim(fnm),false);}
 }
 ff.loadRules(INP.defs_file);
 if(INP.extra_rules!="")
 {
   extra_data+="\tExtra Rules: "+INP.extra_rules.substr(1)+"\n";
   for(const std::string fnm : stringfx::split(INP.extra_rules.substr(1),',')) {std::cout <<"Loading: "<< fnm << "\n"; ff.loadRules(stringfx::trim(fnm),false);}
 }
 if(extra_data!="")
 {
   cout << "Extra FF data:\n";
   cout << extra_data << "\n";
 }
 cin >> tempv;

 IncompleteForceField protff(INP.protein_ff); protff.loadResidues(INP.resid_ff);

 //ML part (start)
 atom_featurizer=new torchmol::FFOneHotFeaturizer(ff,166);
 std::vector<torch::Tensor> params_dec,params_conv;

 #ifdef SERIAL
  SimpleModel sample_dec(2*atom_featurizer->getFeatureSize(),128);
  //RNNModel sample_dec(2*atom_featurizer->getFeatureSize(),128); //Hidden size=128
  //mytorch::modules::GraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize()/*,F::tanh*/);
  mytorch::modules::IterativeGraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),2,torch::tanh);

  molenvmt = new torchmol::GrowingDeNovoMolecule(nullptr,MAX_MOLSIZE,*atom_featurizer,0);
  molenvmt->setMoleculeFeaturizer(*graphconv);
  dnvreinforce::Policy localpolicy(*sample_dec);
  mypolicy = new dnvreinforce::Reinforce(localpolicy);
  torch::load(sample_dec,PTPATH+"/decsave.pt"); //TODO: Fix path
  torch::load(graphconv,PTPATH+"/convsave.pt"); //TODO: Fix path
  graphconv->eval();
  sample_dec->eval();

  mypolicy->replaceOptimizer(new torch::optim::Adam(mytorch::combineParameters(graphconv->parameters(),sample_dec->parameters()),torch::optim::AdamOptions(1e-3))); //Learning rate
  #else
  std::vector<SimpleModel> sample_decs;
  //std::vector<mytorch::modules::GraphConvolution> graphconvs;
  std::vector<mytorch::modules::IterativeGraphConvolution> graphconvs;

  for(int i=0;i<MAX_THREADS_OMP;i++)
  {
    sample_decs.push_back(SimpleModel(2*atom_featurizer->getFeatureSize(),128));
    //sample_decs.push_back(RNNModel(2*atom_featurizer->getFeatureSize(),128));
    //graphconvs.push_back(mytorch::modules::GraphConvolution(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize()/*,F::tanh*/));
    graphconvs.push_back(mytorch::modules::IterativeGraphConvolution(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),2,torch::tanh));
    //graphconvs.push_back(mytorch::modules::IterativeGraphConvolution(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),2,torch::tanh));

    auto* obj=new torchmol::GrowingDeNovoMolecule(nullptr,MAX_MOLSIZE,*atom_featurizer,0,i);
    molenvmt.push_back(obj);
    obj->setMoleculeFeaturizer(*graphconvs[i]);
  }

  std::vector<dnvreinforce::Policy> localpolicies;
  for(int i=0;i<MAX_THREADS_OMP;i++) localpolicies.push_back(dnvreinforce::Policy(*sample_decs[i]));
  mypolicy.init(localpolicies);
  for(int i=0;i<MAX_THREADS_OMP;i++)
  {
      torch::load(sample_decs[i],PTPATH+"/decsave.pt"); //TODO: Fix path
      torch::load(graphconvs[i],PTPATH+"/convsave.pt"); //TODO: Fix path
      sample_decs[i]->eval();
      graphconvs[i]->eval();
  }

  for(int i=0;i<MAX_THREADS_OMP;i++) mypolicy[i]->replaceOptimizer(new torch::optim::Adam(mytorch::combineParameters(graphconvs[i]->parameters(),sample_decs[i]->parameters()),torch::optim::AdamOptions(1e-3)));
  #endif
 //ML part (end)

 Protein p(INP.fnam,protff,hsres,INP.format);
 cout << "Protein loaded\n";
 std::vector<Atom*> seeds;
 for(std::string s : seednames) seeds.push_back(ff.getAtom(stringfx::trim(s)));


 Molecule* templatemol=(INP.prune)?new Molecule(INP.prunemolfile,ff,"pdb",true):nullptr;
 for(int J : INP.size)
 {
    int K=0;
    std::vector<System> mols=p.generateLigands(ff,INP.num,J,INP.progname+"_"+to_string(J)+".log",seeds,INP.seednum,false,INP.spread,1,INP.osc,300,50,298,"seeds_"+INP.progname+"_"+to_string(J)+".pdb",(INP.flagcodes[0]!='0'),(INP.flagcodes[1]!='0'),(INP.flagcodes[2]!='0'),(INP.flagcodes[3]!='0'),region,templatemol,INP.prune,INP.prunefrac);
    for(System& s : mols)
    {
      Molecule* m=s.getDrugMolecule();
      for(Atom* a : m->getAtoms()) cout <<a->toString()<<"	"<<ff.isSatisfied(a,m->getBondedAtoms(a))<<","<<a->isCyclized()<<"\n";
      std::string effname=(stringfx::indexOf("\%i",INP.progname)!=std::string::npos)?stringfx::replace(INP.progname,"\%i",to_string(J)):INP.progname+"_"+to_string(J);
      m->dumpMol("result_prot_"+effname+"_"+to_string(K)+".pdb",&ff);
      m->dumpMolAsMol2("result_prot_"+effname+"_"+to_string(K)+".mol2",&ff,false,true);
      K++;
    }
 }
}
