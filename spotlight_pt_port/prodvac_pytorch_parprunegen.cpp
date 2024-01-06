#define SERIAL 1 //Optional serial generation
#define STATICDATA 1
//#define NOZMQ 1
#define NOJSON 1
#define USE_TORCH 1

#include "graph/Molecule.hpp"

static std::string PTPATH="../";

#define TRAIN_MODE false //DO NOT CHANGE
//#define SEGMENT 4000
#ifndef SERIAL
  static const int MAX_THREADS_OMP=omp_get_max_threads(); //Change here if needed
#endif

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
static double running_reward=-0.3d;
static const bool LOAD_IN_TRAIN_MODE=true;

static Atom* internalGetPolicyAtom(const ForceField& ff,const std::vector<std::string>& atlst,Atom* src,const std::vector<Atom*>& nhd,int devid=0)
{
#ifndef SERIAL
  torch::NoGradGuard tngg;
#endif

  std::vector<Atom*> chs=ff.toAtomList(atlst);
  torch::Tensor allowed_atomfeats=featurizeAtomTypes(*atom_featurizer,chs);
  #ifdef SERIAL
  torch::Tensor sing=molenvmt->getMoleculeFeatures()[molenvmt->getMolecule()->indexOf(src)].repeat({allowed_atomfeats.size(0),1}); //atom_featurizer->featurize(src,nhd).repeat({allowed_atomfeats.size(0),1});
  #else
  torch::Tensor sing=molenvmt[devid]->getMoleculeFeatures()[molenvmt[devid]->getMolecule()->indexOf(src)].repeat({allowed_atomfeats.size(0),1}); //atom_featurizer->featurize(src,nhd).repeat({allowed_atomfeats.size(0),1});
  #endif
  allowed_atomfeats = torch::cat({allowed_atomfeats,sing},1);
  int retidx = 0;
  if(randomchance>0 && RNG::throwarandompoint(0,1)<randomchance) retidx=(int)(RNG::throwarandompoint(0,1)*chs.size());
  #ifdef SERIAL
  else retidx=mypolicy->select_action(allowed_atomfeats);
  mypolicy->save_reward(0.0f);
  #else
  else retidx=mypolicy[devid]->select_action(allowed_atomfeats);
  mypolicy[devid]->save_reward(0.0f);
  #endif
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
static float VARIANCE_REG_CONST=0.01f,REG_SCALE=12.5f,VARIANCE_REG_SCALE=2.000f;
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
static void endState(Molecule* m,bool succ=false,int devid=0)
{
  //torch::Tensor lastvec=builtmols[builtmols.size()-1];
  //float givenReward=succ?ExtendedSYBAScore(m):FAIL_REWARD;
  float givenReward=succ?SYBAScore(m):(FAIL_REWARD);
  if(succ) cout << "Given molecule "<<builtmols.size()<<" a reward of "<<givenReward<<"\n";
  double disc_rew;
  #ifdef SERIAL
  mypolicy->modifyReward(givenReward);
  disc_rew=mypolicy->finish_episode(/*train=*/TRAIN_MODE);
  #else
  mypolicy[devid]->modifyReward(givenReward);
  disc_rew=mypolicy[devid]->finish_episode(/*train=*/false);
  #endif

  if(succ)
  {
#ifdef SERIAL
  	running_reward=running_reward*(1.0-1.0/NUM_AVG)+givenReward/NUM_AVG;
#endif
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
		zmqio::initThreads(cout,MAX_THREADS_OMP);
	#else
		zmqio::initThreads(cout);
	#endif
}

int main(int argc,char** argv)
{
#ifndef SERIAL
  torch::NoGradGuard tngg;
#endif
  USE_PHYSICS=1;
  FORCE_CHARGES=1;
  FFCHARGES=true;
  QUIET=1;
  SYNTHGET=1;
  SYNTHPENALTY=350;
  WEIGHTED_ATOM_SELECTOR_FUNCTION=internalGetPolicyAtom;
  //IPADDR="127.0.0.1";
  //PORT=5555;

  constexpr int HEADS=8;
  constexpr bool concatGAT=false;

  constexpr int outsize=(concatGAT)?(HEADS+1):2;
  
  setupNetwork();
  RNG::init();
#ifndef SERIAL
  std::cout << "OMP Thread Count: "<<MAX_THREADS_OMP<<"\n";
  std::cout.flush();
#endif
  
  //const int size=(argc>1)?std::stoi(argv[1]):80,nmol=(argc>2)?std::stoi(argv[2]):250,vari=(argc>3)?std::stoi(argv[3]):0;
  const std::string molfile=std::string(argv[1]);
  const float cutfrac=std::stof(argv[2]);
  const int nmol=std::stoi(argv[3]),vari=std::stoi(argv[4]);
  const bool KEEPHYDS=true;
  ff.loadBondParameters("/home/venkata/dnv/data/itps/bondtypes.itp","/home/venkata/dnv/data/itps/angletypes.itp","/home/venkata/dnv/data/itps/dihedraltypes.itp","/home/venkata/dnv/data/itps/impropertypes.itp");
  ff.loadRules("/home/venkata/dnv/data/definitions.data");

  Molecule* loadmol=new Molecule(molfile,ff,"pdb",true);
  loadmol->type=0;

  atom_featurizer=new torchmol::FFOneHotFeaturizer(ff);

  #ifdef SERIAL
  SimpleModel sample_dec(outsize*atom_featurizer->getFeatureSize(),128);
  //RNNModel sample_dec(2*atom_featurizer->getFeatureSize(),128); //Hidden size=128

  mytorch::modules::IterativeGraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),4,torch::tanh);
  //mytorch::modules::GraphConvolution graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize()/*,F::tanh*/);
  //mytorch::modules::GraphAttention graphconv(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),HEADS,torch::relu/*,torch::tanh*/,concatGAT);


  molenvmt = new torchmol::GrowingDeNovoMolecule(nullptr,MAX_MOLSIZE,*atom_featurizer,0);
  molenvmt->setMoleculeFeaturizer(*graphconv);
  dnvreinforce::Policy localpolicy(*sample_dec);
  mypolicy = new dnvreinforce::Reinforce(localpolicy);

  #ifdef SEGMENT
  	torch::load(sample_dec,PTPATH+"/decsave_"+std::to_string(SEGMENT)+".pt");
  	torch::load(graphconv,PTPATH+"/convsave_"+std::to_string(SEGMENT)+".pt");
  #else
  	torch::load(sample_dec,PTPATH+"/decsave.pt");
  	torch::load(graphconv,PTPATH+"/convsave.pt");
  #endif

  graphconv->eval();
  sample_dec->eval();

  mypolicy->replaceOptimizer(new torch::optim::Adam(mytorch::combineParameters(graphconv->parameters(),sample_dec->parameters()),torch::optim::AdamOptions(1e-3))); //Learning rate

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
  #else
  std::vector<SimpleModel> sample_decs;
  //std::vector<RNNModel> sample_decs;
  //std::vector<mytorch::modules::GraphConvolution> graphconvs;
  std::vector<mytorch::modules::IterativeGraphConvolution> graphconvs;
  //std::vector<mytorch::modules::GraphAttention> graphconvs;

  for(int i=0;i<MAX_THREADS_OMP;i++)
  {
    sample_decs.push_back(SimpleModel(2*atom_featurizer->getFeatureSize(),128));
    //sample_decs.push_back(RNNModel(2*atom_featurizer->getFeatureSize(),128));

    //graphconvs.push_back(mytorch::modules::GraphConvolution(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize()/*,F::tanh*/));
    graphconvs.push_back(mytorch::modules::IterativeGraphConvolution(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),4,torch::tanh));
    //graphconvs.push_back(mytorch::modules::GraphAttention(atom_featurizer->getFeatureSize(),atom_featurizer->getFeatureSize(),HEADS,torch::relu/*,torch::tanh*/,concatGAT));

    auto* obj=new torchmol::GrowingDeNovoMolecule(nullptr,MAX_MOLSIZE,*atom_featurizer,0,i);
    molenvmt.push_back(obj);
    obj->setMoleculeFeaturizer(*graphconvs[i]);
  }

  std::vector<dnvreinforce::Policy> localpolicies;
  for(int i=0;i<MAX_THREADS_OMP;i++) localpolicies.push_back(dnvreinforce::Policy(*sample_decs[i]));
  mypolicy.init(localpolicies);
  for(int i=0;i<MAX_THREADS_OMP;i++)
  {
	#ifdef SEGMENT
  	torch::load(sample_decs[i],PTPATH+"/decsave_"+std::to_string(SEGMENT)+".pt");
  	torch::load(graphconvs[i],PTPATH+"/convsave_"+std::to_string(SEGMENT)+".pt");
	#else
  	torch::load(sample_decs[i],PTPATH+"/decsave.pt");
  	torch::load(graphconvs[i],PTPATH+"/convsave.pt");
	#endif
      sample_decs[i]->eval();
      graphconvs[i]->eval();
  }

  for(int i=0;i<MAX_THREADS_OMP;i++) mypolicy[i]->replaceOptimizer(new torch::optim::Adam(mytorch::combineParameters(graphconvs[i]->parameters(),sample_decs[i]->parameters()),torch::optim::AdamOptions(1e-3)));
  #endif

  std::ofstream of; of.open("dnv_prod.log",ios::out);
  #ifdef SERIAL
  int tid=0;
  for(int i=0;i<nmol;i++)
  #else
  #pragma omp parallel
  for(int i=omp_get_thread_num();i<nmol;i+=MAX_THREADS_OMP)
  #endif
  {
    int tid=omp_get_thread_num();
    System sys(ff);

    int esize=loadmol->getEffectiveSize()+((vari)?(int)throwarandompoint(-vari-1,vari+1):0);
    //Atom* root=ff.getRandomAtomByType();
    Molecule* m =new Molecule(loadmol); m->type=0;
    if(!m->makeRandomPruneFragment(&ff,cutfrac,KEEPHYDS))
    {
      cerr << "ERR: Could not make a prune fragment!!!\n";
      assert(1!=1);
    }
    else
    {
      std::cout << "Free atoms: "<<m->freeatoms.size()<<"\n";
      for(Atom* fa : m->freeatoms) std::cout<<"\t" << fa->toString() << "\n";
      m->dumpMol("created_frag_"+std::to_string(tid)+".pdb",&ff);
      m->describeStructure();
    }
    sys.addMolecule(m);

#ifndef SERIAL
    molenvmt[tid]->replaceMolecule(m);
    try {m->samplegrow(of,sys,esize,40,100,20,false,false,false,false,NULLMOL,0,tid);}
#else
    molenvmt->replaceMolecule(m);
    try {m->samplegrow(of,sys,esize,40,100,20);}
#endif
    catch(NoBondsAvailableException ex) 
    {
#ifdef SERIAL
      i--;
      endState(m,false);
#else
      i-=MAX_THREADS_OMP;
      endState(m,false,tid);
#endif
      continue;
    }
    if(m->getEffectiveSize()<esize || !ff.isSatisfied(m->getAtoms()[0],m->getBondedAtoms(m->getAtoms()[0])))
    {
      cout << "Failed!\n";

    #ifdef SERIAL
      i--;
      endState(m,false);
     #else
      i-=MAX_THREADS_OMP;
      endState(m,false,tid);
     #endif
      continue;
    }
	m->calculateBondOrders2(ff); m->pE=0.0;
	netcomm::processRequest(of,m,tid);
	if(!m)
    {
      cout << "Bad molecule by SMILES\n";
     #ifdef SERIAL
      i--;
      endState(m,false);
     #else
      i-=MAX_THREADS_OMP;
      endState(m,false,tid);
     #endif
      continue;
    }
    //builtmols.push_back(getAtomTypeWeights(featurizeAtomTypes(*atom_featurizer,m->getAtoms())));

  #ifdef SERIAL
    endState(m,true);
  #else
    endState(m,true,tid);
  #endif

    cout << "Done training\n";
    m->dumpMol("result_prodvac_"+to_string(i)+".pdb",&ff);
#ifdef SERIAL
    if(TRAIN_MODE && i%50==0)
    {
    	torch::save(sample_dec,PTPATH+"/decsave.pt");
    	torch::save(graphconv,PTPATH+"/convsave.pt");

    	torch::save(sample_dec,PTPATH+"/decsave_"+to_string(i)+".pt");
    	torch::save(graphconv,PTPATH+"/convsave_"+to_string(i)+".pt");

    	cout << "Saved model till now\n";
    }
#endif
    std::cout << "Running Reward: " << running_reward << "\n";
    if(running_reward>TARGET_REWARD && i>NUM_AVG && TRAIN_MODE) break;
    //delete m; // Optional
  }
  #ifndef SERIAL
  #pragma omp barrier
  #else
  torch::save(sample_dec->parameters(),"../test_dec.pt");
  torch::save(graphconv->parameters(),"../test_conv.pt");
  std::ofstream outf; outf.open("tensordump.dat");
  outf << "Decoder parameters\n\n";
  for(auto& parm : sample_dec->parameters()) outf << parm << "\n";
  outf << "Convolver parameters\n\n";
  for(auto& parm : graphconv->parameters()) outf << parm << "\n";
  outf.close();
  //Using torch save
  torch::save(sample_dec,PTPATH+"/decsave.pt");
  torch::save(graphconv,PTPATH+"/convsave.pt");
  cout << "Saved final model\n";
  #endif
  cout << "Done with "<<nmol<<" molecules.\n";
  of.close();
}
