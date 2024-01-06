#define SERIAL 1
#define STATICDATA 1
//#define NOZMQ 1
#define NOJSON 1
#define USE_TORCH 1

#include "graph/Molecule.hpp"

static std::string PTPATH="../";

//set to false to stop - TRAIN MODE DISABLES LOADING OF PREVIOUS PARAMETERS
#define TRAIN_MODE true
//#define SEGMENT 1900

// A fixed size queue-like list that automatically maintains size by removing oldest elements
template<class T> class FixedSizeVector
{
	T* holding;
	int sizelim=0,trueSize=0;
	int ind=0;
	bool filled=false;
public:
	FixedSizeVector(const int sz)
	{
		sizelim=sz;
		holding=new T[sz];
		ind=0;
	}
	~FixedSizeVector() {delete[] holding;}
	
	/*inline virtual int correctedIndex(int idx) const {return idx%sizelim;}
	inline T& operator[](int idx) {return holding[correctedIndex(idx)];}
	inline const T& operator[](int idx) const {return holding[correctedIndex(idx)];}*/
	
	void push_back(const T& obj)
	{
		holding[ind++]=obj;
		trueSize++;
		if(ind==sizelim) {ind=0; filled=true;}
	}
	inline std::vector<T> toStandardVector() const {return std::vector<T>(holding,holding+this->size());}
	inline operator std::vector<T>() const {return toStandardVector();}
	inline torch::Tensor getLastTensor() const {return holding[ind?(ind-1):(sizelim-1)];}
	inline int size() const {return ((filled)?sizelim:ind);}
	inline int actual_size() const {return trueSize;}
};

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

static const int BLOCK_SIZE=50,NUM_AVG=50;
static std::vector<torch::Tensor> builtmols;
static const double TARGET_REWARD=0.20;


static ForceField ff("/home/venkata/dnv/data/final_ff_parameters.ffin","/home/venkata/dnv/data/categories.data");
static float randomchance=0.0;
static double running_reward=-0.45d;
static const bool LOAD_IN_TRAIN_MODE=true;

static Atom* internalGetPolicyAtom(const ForceField& ff,const std::vector<std::string>& atlst,Atom* src,const std::vector<Atom*>& nhd,int devid=0)
{
  std::vector<Atom*> chs=ff.toAtomList(atlst);
  torch::Tensor allowed_atomfeats=featurizeAtomTypes(*atom_featurizer,chs);
  torch::Tensor sing=molenvmt->getMoleculeFeatures()[molenvmt->getMolecule()->indexOf(src)].repeat({allowed_atomfeats.size(0),1}); //atom_featurizer->featurize(src,nhd).repeat({allowed_atomfeats.size(0),1});
  allowed_atomfeats = torch::cat({allowed_atomfeats,sing},1);
  int retidx = 0;
  if(randomchance>0 && RNG::throwarandompoint(0,1)<randomchance) retidx=(int)(RNG::throwarandompoint(0,1)*chs.size());
  else retidx=mypolicy->select_action(allowed_atomfeats);
  mypolicy->save_reward(0.0f);
  return chs[retidx];
}


static float fractionCyclizedScore(Molecule* mol)
{
  int nCyc=0;
  for(Atom* at : mol->getAtoms())
  {
    if(at->isCyclized()) nCyc++;
  }
  return (float)nCyc/(float)mol->getSize();
}
inline static float SYBAScore(Molecule* mol) {return (mol->pE>150)?-1.0f:1.0f;}

// Regularization scoring START
static int MOL_REG_CONST=125,MOL_REG_THRESH=5;
static float VARIANCE_REG_CONST=0.01f,REG_SCALE=12.5f,VARIANCE_REG_SCALE=2.250f;
static torch::Tensor getAtomTypeWeights(torch::Tensor t)
{
	torch::Tensor cnts=torch::sum(t,0);
	float norm=torch::sum(t).item<float>();
	return cnts/norm;
}
static float getVarianceScore(std::vector<torch::Tensor> tensesin,int take=BLOCK_SIZE) // tenses is the set of feature-vectors (here a feature-vector is just the normalized atom-type count) for each molecule accepted till now
{
	take=::min(take,(int)tensesin.size()-1);
	if(tensesin.size()<MOL_REG_THRESH) return 0.0f;
	if(!take) take=tensesin.size()-1;
	
	cerr << tensesin.size() << " as input size and take = "<<take<<"\n";
	int ubound=(tensesin.size()-2),lbound=tensesin.size()-2-take; //'take' has been corrected to not exceed len(list)-1
	cerr << "Upper bound: "<< ubound <<" and Lower bound: "<< lbound << "\n";
	
	std::vector<torch::Tensor> tenses;
	for(int idx1=ubound;idx1>lbound;idx1--) {tenses.push_back(tensesin[idx1]); cerr << idx1 <<", ";}
	cerr <<". So "<< tenses.size() << " is stack size.\n";
	
	torch::Tensor tsstack=torch::stack(tenses);
	cout << mytorch::getShape(tsstack) << " ";
	torch::Tensor mytens=tenses[tenses.size()-1];
	/*
		// Use cosine similarity to mean feature vector
		
		torch::Tensor stackmean=torch::mean(tsstack,0); // Has Shape (166,) - Is correct (AFAIK)
		
		float dotprod=torch::dot(stackmean,mytens).item<float>();
		float norms=torch::norm(stackmean).item<float>()*torch::norm(mytens).item<float>();
		return -(dotprod/norms)*VARIANCE_REG_SCALE;
	*/
	
	//Use L2 norm of pairwise cosine similarities
	torch::Tensor norms=torch::linalg::vector_norm(tsstack,/*ord=*/2,/*dim=*/1,/*keepdim=*/false,torch::kFloat32);
	torch::Tensor dotvs=torch::mm(tsstack,mytens.unsqueeze(1)).squeeze()/torch::norm(mytens);
	std::cout << mytorch::getShape(norms) << " as norm shape\n";
	dotvs=dotvs/norms;
	std::cout << "\n\n" << dotvs << " as dot prods\n";
	return -::sqrt(torch::mean(dotvs*dotvs).item<float>())*VARIANCE_REG_SCALE;
}
// Regularization scoring END

static float ExtendedSYBAScore(Molecule* mol) {return (mol->pE>150)?-1.0f:(1.0f+getVarianceScore(builtmols));}
static void endState(Molecule* m,bool succ=false)
{
  //torch::Tensor lastvec=builtmols[builtmols.size()-1];
  //float givenReward=succ?ExtendedSYBAScore(m):FAIL_REWARD;
  float givenReward=succ?ExtendedSYBAScore(m):(FAIL_REWARD);
  if(succ) cout << "Given molecule "<<builtmols.size()<<" a reward of "<<givenReward<<"\n";
  mypolicy->modifyReward(givenReward);
  double disc_rew=mypolicy->finish_episode(/*train=*/TRAIN_MODE);
  if(succ)
  {
  	running_reward=running_reward*(1.0-1.0/NUM_AVG)+givenReward/NUM_AVG;
	std::cout << "Given Reward: "<<givenReward<<"\n";
	std::cout << "Discounted Reward: "<<disc_rew<<"\n";
  }
  if(!succ && m) delete m;
}

static void setupNetwork()
{
	int K[2]; K[0]=0; K[1]=0;
	netcomm::genff=&ff;
	netcomm::badsmileslog.open("badsmiles.log");
	netcomm::failures=&K[0];
	
	IPADDR="127.0.0.1";
	PORT=5555;
	
	#ifndef SERIAL
		zmqio::initThreads(cout,omp_get_max_threads());
	#else
		zmqio::initThreads(cout);
	#endif
}

int main(int argc,char** argv)
{
  USE_PHYSICS=1;
  FORCE_CHARGES=1;
  FFCHARGES=true;
  QUIET=1;
  SYNTHGET=1;
  SYNTHPENALTY=350;
  WEIGHTED_ATOM_SELECTOR_FUNCTION=internalGetPolicyAtom;
  //IPADDR="127.0.0.1";
  //PORT=5555;

  constexpr int HEADS=6;
  constexpr bool concatGAT=false;

  constexpr int outsize=(concatGAT)?(HEADS+1):2;
  
  setupNetwork();
  RNG::init();
  
  const int size=(argc>1)?std::stoi(argv[1]):80,nmol=(argc>2)?std::stoi(argv[2]):250,vari=(argc>3)?std::stoi(argv[3]):0;
  ff.loadBondParameters("/home/venkata/dnv/data/itps/bondtypes.itp","/home/venkata/dnv/data/itps/angletypes.itp","/home/venkata/dnv/data/itps/dihedraltypes.itp","/home/venkata/dnv/data/itps/impropertypes.itp");
  ff.loadRules("/home/venkata/dnv/data/definitions.data");

  atom_featurizer=new torchmol::FFOneHotFeaturizer(ff);
  SimpleModel sample_dec(outsize*atom_featurizer->getFeatureSize(),128);
  //RNNModel sample_dec(2*atom_featurizer->getFeatureSize(),128); //Hidden size=128
  
  //mytorch::modules::IterativeGraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),2,torch::tanh);
  //mytorch::modules::GraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize()/*,F::tanh*/);
  mytorch::modules::GraphAttention graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),HEADS,torch::relu/*,torch::tanh*/,concatGAT);

  molenvmt = new torchmol::GrowingDeNovoMolecule(nullptr,MAX_MOLSIZE,*atom_featurizer,0);
  molenvmt->setMoleculeFeaturizer(*graphconv);

  std::vector<torch::Tensor> params_dec,params_conv;
  dnvreinforce::Policy localpolicy(*sample_dec);
  mypolicy = new dnvreinforce::Reinforce(localpolicy);
  
 
  /*torch::load(params_dec,"../test_dec.pt"); torch::load(params_conv,"../test_conv.pt");
  cout << "INFO: "<< params_dec.size() << " and " << params_conv.size() << " as loaded parameter counts for decoder and convolver respectively.\n";
  int myctr=0;
  for(auto& parm : sample_dec->parameters()) parm = params_dec[myctr++];
  myctr=0;
  for(auto& parm : graphconv.parameters()) parm = params_conv[myctr++];
  cout << "Model weights loaded\n";*/
  
  if(!TRAIN_MODE || LOAD_IN_TRAIN_MODE)
  {
	#ifdef SEGMENT
  	torch::load(sample_dec,PTPATH+"/decsave_"+std::to_string(SEGMENT)+".pt");
  	torch::load(graphconv,PTPATH+"/convsave_"+std::to_string(SEGMENT)+".pt");
	#else
  	torch::load(sample_dec,PTPATH+"/decsave.pt");
  	torch::load(graphconv,PTPATH+"/convsave.pt");
	#endif
  }
  
  
  mypolicy->replaceOptimizer(new torch::optim::Adam(mytorch::combineParameters(sample_dec->parameters(),graphconv->parameters()),torch::optim::AdamOptions(0.95e-3))); //Learning rate

  //torchmol::weightedAtomSelector=new torchmol::AtomSelector(ff,*atom_featurizer,sample_dec);

  std::vector<Atom*> seeds;
  seeds.push_back(ff.getAtom("CT1"));
  seeds.push_back(ff.getAtom("CT2"));
  seeds.push_back(ff.getAtom("CA"));
  seeds.push_back(ff.getAtom("CA"));
  seeds.push_back(ff.getAtom("NX"));
  seeds.push_back(ff.getAtom("CAS"));
  seeds.push_back(ff.getAtom("NX"));
  seeds.push_back(ff.getAtom("CPT"));
  seeds.push_back(ff.getAtom("CPT"));
  seeds.push_back(ff.getAtom("NY"));
  seeds.push_back(ff.getAtom("CY"));
  seeds.push_back(ff.getAtom("CT2x"));
  seeds.push_back(ff.getAtom("CT2x"));
  seeds.push_back(ff.getAtom("CT1x"));
  seeds.push_back(ff.getAtom("CT1x"));
  seeds.push_back(ff.getAtom("CP2"));
  seeds.push_back(ff.getAtom("CE1A"));
  seeds.push_back(ff.getAtom("CE2A"));
  seeds.push_back(ff.getAtom("PL"));
  seeds.push_back(ff.getAtom("SL"));
  seeds.push_back(ff.getAtom("C3"));
  seeds.push_back(ff.getAtom("OH1"));
  seeds.push_back(ff.getAtom("OC"));

  seeds[2]->mustCycl=true;
  seeds[3]->mustCycl=true;
  seeds[4]->mustCycl=true;
  seeds[5]->mustCycl=true;
  seeds[6]->mustCycl=true;
  seeds[7]->mustCycl=true;
  seeds[8]->mustCycl=true;
  seeds[9]->mustCycl=true;
  seeds[10]->mustCycl=true;
  seeds[11]->mustCycl=true;
  seeds[12]->mustCycl=true;
  seeds[13]->mustCycl=true;
  seeds[14]->mustCycl=true;
  seeds[15]->mustCycl=true;
  seeds[20]->mustCycl=true;

  //seeds.push_back(ff.getAtom("C6x"));
  //seeds[0]->mustCycl=true;

  std::ofstream of; of.open("dnv_prod.log",ios::out);
  for(int i=0;i<nmol;i++)
  {
    System sys(ff);
    //Atom* root=ff.getRandomNonHAtomByType();
    Atom* root=new Atom(*randomSelect(seeds));
    cout << "Picked seed\n";
    cout << root->toString() << "as seed\n";
    //Atom* root=ff.getRandomAtomByType();
    Molecule* m =new Molecule(root);
    sys.addMolecule(m);
    molenvmt->replaceMolecule(m);
    int esize=size+((vari)?(int)throwarandompoint(-vari-1,vari+1):0);
    try {m->samplegrow(of,sys,esize,40,100,20);}
    catch(NoBondsAvailableException ex) 
    {
      i--; 
      endState(m,false);
      continue;
    }
    if(m->getEffectiveSize()<esize || !ff.isSatisfied(m->getAtoms()[0],m->getBondedAtoms(m->getAtoms()[0]))) {cout << root->toString()<<"\t"; cout << "Failed!"; i--; endState(m,false); continue;}
	m->calculateBondOrders2(ff); m->pE=0.0;
	netcomm::processRequest(of,m,0);
	if(!m) {cout << "Bad molecule by SMILES\n"; i--; endState(m,false); continue;}
    cout << root->toString()<<"\tSuccess!\n";
    builtmols.push_back(getAtomTypeWeights(featurizeAtomTypes(*atom_featurizer,m->getAtoms())).detach());
    endState(m,true);
    cout << "Done training\n";
    m->dumpMol("result_prodvac_"+to_string(i)+".pdb",&ff);
    if(TRAIN_MODE && i%50==0)
    {
    	torch::save(sample_dec->parameters(),"../test_dec.pt");
    	torch::save(graphconv->parameters(),"../test_conv.pt");
    	torch::save(sample_dec,PTPATH+"/decsave.pt");
    	torch::save(graphconv,PTPATH+"/convsave.pt");

    	torch::save(sample_dec->parameters(),"../test_dec_"+to_string(i)+".pt");
    	torch::save(graphconv->parameters(),"../test_conv_"+to_string(i)+".pt");
    	torch::save(sample_dec,PTPATH+"/decsave_"+to_string(i)+".pt");
    	torch::save(graphconv,PTPATH+"/convsave_"+to_string(i)+".pt");

    	cout << "Saved model till now\n";
    }
    std::cout << "Running Reward: " << running_reward << "\n";
    if(running_reward>TARGET_REWARD && i>NUM_AVG && TRAIN_MODE) break;
    //delete m; // Optional
  }
  cout << "Done with "<<nmol<<" molecules.\n";
  of.close();
  torch::save(sample_dec->parameters(),"../test_dec.pt");
  torch::save(graphconv->parameters(),"../test_conv.pt");
  
  std::ofstream outf; outf.open("tensordump.dat");
  outf << "Decoder parameters\n\n";
  for(auto& parm : sample_dec->parameters()) outf << parm << "\n";
  outf << "Convolver parameters\n\n";
  for(auto& parm : graphconv->parameters()) outf << parm << "\n";
  outf.close();
  //((const mytorch::TorchscriptModule&)(mypolicy->getPolicy().getModelModule())).getTorchscriptModule().save("../test.ts");
  
  //Using torch save
  torch::save(sample_dec,PTPATH+"/decsave.pt");
  torch::save(graphconv,PTPATH+"/convsave.pt");
  cout << "Saved final model\n";
}
