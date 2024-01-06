#ifndef INCLUDED_TORCHMOL_DNV
#define INCLUDED_TORCHMOL_DNV 1
#include "support/mytorch.h"
//#include "graph/Molecule.hpp"
#include <iostream>
// using std::cout;
typedef int64_t TorchSizeElement;
using namespace torch::data::datasets;
using namespace torch::indexing;
//typedef at::NoGradGuard torch::NoGradGuard;
namespace F = torch::nn::functional;

typedef std::pair<torch::Tensor,torch::Tensor> TorchGraph;
namespace torchmol
{
  class AtomFeaturizer : public torch::nn::Module
  {
  protected:
    int feat_dim=-1;
  public:
    AtomFeaturizer(int fdim) : torch::nn::Module() {feat_dim=fdim;}

    virtual torch::Tensor featurize(Atom*,const std::vector<Atom*>&) const =0;
    inline torch::Tensor /*operator()*/forward(Atom* at,const std::vector<Atom*>& bnd) const {return this->featurize(at,bnd);}
    inline torch::Tensor /*operator()*/forward(Atom* at,Molecule* m) const {return this->featurize(at,m->getBondedAtoms(at));}
    inline int getFeatureSize() const {return feat_dim;}
  };
  class FFOneHotFeaturizer : public AtomFeaturizer
  {
    const ForceField* myff;
    std::map<std::string,int> atom_inds;
    std::map<int,std::string> atom_names;
  public:
    FFOneHotFeaturizer(const ForceField* ff,int force_dim=-1) : AtomFeaturizer((force_dim<0)?ff->getAtomTypes().size():force_dim) {myff=ff; this->refreshIndices();}
    FFOneHotFeaturizer(const ForceField& ff,int force_dim=-1) : FFOneHotFeaturizer(&ff,force_dim) {}

    void refreshIndices()
    {
      atom_inds=std::map<std::string,int>();
      atom_names=std::map<int,std::string>();
      const auto& myatoms = myff->getAtomTypes();
      for(int i=0;i<myatoms.size();i++)
      {
        atom_inds.emplace(myatoms[i]->toString(),i);
        atom_names.emplace(i,myatoms[i]->toString());
      }
    }

    torch::Tensor featurize(Atom* cent,const std::vector<Atom*>& nhd=std::vector<Atom*>()) const override
    {
      /*std::vector<torch::Tensor> onehots;
      for(int i=0;i<m->getSize();i++)
      {
        Atom* atm = m->getAtom(i);
        torch::Tensor encd=torch::zeros(feat_dim);
        encd[atom_inds.at(atm->toString())]=1.;
        onehots.push_back(encd);
      }
      torch::Tensor ret = torch::zeros({MAX_MOLSIZE,feat_dim});
      ret.index_put({Slice(0,onehots.size()),Slice(0,feat_dim)},torch::stack(onehots));
      return ret;*/
      torch::Tensor ret=torch::zeros(feat_dim);
      int fidx=atom_inds.at(cent->toString());
      if(fidx<feat_dim) ret[fidx]=1.;
      return ret;
    }
  };

  class NoMoleculeLevelFeaturizerException : public std::exception {};

  class GrowingDeNovoMolecule
  {
    Molecule* mol;
    int maxsize=0;
    int focus=-1;
    const AtomFeaturizer* feat;

    torch::jit::script::Module mod;
    mytorch::TorchModuleWrapper* molfeat=nullptr;
    bool molfeatactive=false;

    std::vector<torch::Tensor> nodes;
    torch::Tensor adjmat;
    std::vector<torch::Tensor> oldnodes;
    torch::Tensor oldadjmat;
  public:
    const int mydev;

  public:
    GrowingDeNovoMolecule(Molecule* root,int maxs,const AtomFeaturizer& atom_feat,int foc=-1,int pindev=0) : mydev(pindev)
    {
      mol=root;
      maxsize=maxs;
      feat=&atom_feat;
      if(foc<0) foc=mol->getSize()-1;
      focus=foc;
    }

    inline const Molecule* getMolecule() const {return mol;}
    inline Molecule* getMolecule() {return mol;}
    inline void replaceMolecule(Molecule* nmol)
    {
        mol=nmol;
        focus=mol->getSize()-1;
        this->getGraphInput(/*force_recalc=*/true);
    }

    inline const std::vector<torch::Tensor>& getNodes() const {return nodes;}
    inline const torch::Tensor& getAdjacencyMatrix() const {return adjmat;}
    
    TorchGraph getGraphInput(bool force_recalc=false)
    {
      if(nodes.size()!=mol->getSize() || force_recalc)
      {
        nodes=std::vector<torch::Tensor>();
        for(Atom* at : mol->getAtoms()) nodes.push_back(feat->featurize(at,mol->getBondedAtoms(at)));

        //Filling in the AdjMat
        std::vector<int64_t> adjshape; adjshape.push_back(maxsize); adjshape.push_back(maxsize);
        adjmat = torch::zeros(adjshape);
        const auto& mybonds = mol->getBonds();
        for(int i=0;i<mybonds.size();i++)
        {
          for(int j=0;j<mybonds[i].size();j++) adjmat.index_put_({i,mybonds[i][j]},1.);
        }
      }

      torch::Tensor nodefeat = torch::stack(nodes);
      nodefeat = torch::cat({nodefeat,torch::zeros({maxsize-(long)(nodes.size()),mytorch::getShape(nodefeat)[1]})});
      return make_pair(nodefeat,adjmat);
    }

    void updateLastAddedAtom()
    {
      int lind=mol->getSize()-1;
      Atom* lat=mol->getAtom(lind);
      nodes.push_back(feat->featurize(lat,mol->getBondedAtoms(lat)));
      const auto& mybonds = mol->getBonds();
      assert(mol->getSize()<MAX_MOLSIZE);
      for(int i=0;i<mybonds[lind].size();i++)
      {
        adjmat.index_put_({mybonds[lind][i],lind},1.);
        adjmat.index_put_({lind,mybonds[lind][i]},1.);
      }
    }
    
    void startTrials()
    {
    	oldnodes=std::vector<torch::Tensor>(); for(auto& tens : nodes) oldnodes.push_back(tens.clone());
    	oldadjmat=adjmat.clone();
    }
    
    void rollback()
    {
    	nodes=std::vector<torch::Tensor>(); for(auto& tens : oldnodes) nodes.push_back(tens.clone());
    	adjmat=oldadjmat.clone();
    }
    inline void commit() {}

    inline const AtomFeaturizer& getAtomFeaturizer() const {return *feat;}

    inline void setMoleculeFeaturizer(mytorch::TorchModuleWrapper& mfeat) {molfeat=&mfeat; molfeatactive=true;}
    inline bool loadMoleculeFeaturizerFromTorchscript(const std::string& tsfile)
    {
      try
      {
        molfeat = new mytorch::TorchscriptModule(tsfile);
        return (molfeatactive=true);
      }
      catch(const std::exception& e) {return false;}
    }
    inline bool hasMoleculeFeaturizer() const {return molfeatactive;}
    inline const auto& getMoleculeFeaturizer() const {return *molfeat;}
    torch::Tensor getMoleculeFeatures()
    {
      if(!molfeatactive) throw NoMoleculeLevelFeaturizerException();
      //molfeat must take in nodes (hot-encoded or atom-level encoded) and produce more complex features (such as by convolution)
      auto graph = this->getGraphInput();
      
      torch::Tensor nodefeat=get<0>(graph),adjmat=get<1>(graph).clone();
      torch::Tensor encoded_nodes;
      std::vector<torch::Tensor> inps; inps.push_back(nodefeat.unsqueeze(0)); inps.push_back(adjmat.unsqueeze(0));
      //inps[0].set_requires_grad(true); inps[1].set_requires_grad(true);

			encoded_nodes = molfeat->forward(inps);
			return encoded_nodes.squeeze(0);
    }
  };

  static torch::Tensor featurizeAtomTypes(const AtomFeaturizer& feat, const std::vector<Atom*>& ins)
  {
    std::vector<torch::Tensor> ret;
    for(Atom* at : ins) {ret.push_back(feat.featurize(at,std::vector<Atom*>()));}
    return torch::stack(ret);
  }

  /**@brief Small wrapper to use ForceField::listAtomsByRule to get a list of allowed atoms */
  static std::vector<Atom*> getAllowedAtoms(Molecule* loadmol,Atom* selat,const std::vector<Atom*>& selatnhd,const ForceField& ff)
  {
    std::vector<std::string> allowed_atoms = ff.listAtomsByRule(selat,selatnhd);
    /*cout << "Selected Atom: "<<selat->toString()<<"\n";
    for(std::string allat : allowed_atoms) cout << allat<<" ";
    cout << "\n";*/
    return ff.toAtomList(allowed_atoms); //Does NOT copy the Atom* objects
  }

  class AtomWeightGenerator : public torch::nn::Module
  {
    torch::jit::script::Module dec;
    const ForceField* ff;
  public:
    AtomWeightGenerator(const ForceField& myff,torch::jit::script::Module& d) {ff=&myff; dec=d;}

    torch::Tensor forward(GrowingDeNovoMolecule& gmol,int focal) //gmol is guaranteed to have atom-level featurizer
    {
      torch::Tensor encoded_nodes=gmol.getMoleculeFeatures();
			
			Atom* selat = gmol.getMolecule()->getAtom(focal);
			std::vector<Atom*> selatnhd=gmol.getMolecule()->getBondedAtoms(selat);
			std::vector<Atom*> allowed_atoms = getAllowedAtoms(gmol.getMolecule(),selat,selatnhd,*ff);

			torch::Tensor allowed_atomfeats=featurizeAtomTypes(gmol.getAtomFeaturizer(),allowed_atoms);
			torch::Tensor sing=encoded_nodes.select(1,gmol.getMolecule()->indexOf(selat)).repeat({allowed_atomfeats.size(0),1});
			allowed_atomfeats = torch::cat({allowed_atomfeats,sing},1);

      torch::Tensor selwts = mytorch::applyTSModule(dec,allowed_atomfeats).squeeze(1); //sample_dec.forward(encoded_nodes).squeeze();
      return F::softmax(selwts,-1);
    }
  };

  class AtomSelector
  {
    torch::jit::script::Module dec;
    const AtomFeaturizer* atom_feat;
    const ForceField* ff;
  public:
    AtomSelector(const ForceField& myff,const AtomFeaturizer& atfeat, torch::jit::script::Module& d) {ff=&myff; atom_feat=&atfeat; dec=d;}

    Atom* select(Atom* src,const std::vector<Atom*>& bnd,const std::vector<Atom*>& chs) //gmol is guaranteed to have atom-level featurizer
    {
      torch::Tensor allowed_atomfeats=featurizeAtomTypes(*atom_feat,chs);
      torch::Tensor sing=atom_feat->featurize(src,bnd).repeat({allowed_atomfeats.size(0),1});

			//torch::Tensor allowed_atomfeats=featurizeAtomTypes(gmol.getAtomFeaturizer(),allowed_atoms);
			//torch::Tensor sing=encoded_nodes.select(1,gmol.getMolecule()->indexOf(selat)).repeat({allowed_atomfeats.size(0),1});
			allowed_atomfeats = torch::cat({allowed_atomfeats,sing},1);
			cout << mytorch::getShape(allowed_atomfeats) <<" as final feature shape\n";

      torch::Tensor selwts = mytorch::applyTSModule(dec,allowed_atomfeats).squeeze(1); //sample_dec.forward(encoded_nodes).squeeze();
      torch::Tensor finwts = F::softmax(selwts,-1);

      int idx=torch::multinomial(finwts,1,true).squeeze().item<int>();
      cout << finwts << "\n";
      return chs[idx];
    }
  };

  static AtomSelector* weightedAtomSelector=nullptr;
  class NoAtomSelectorLoadedException : public std::exception {};
  static Atom* pickAtomByTorchWeights(const ForceField& ff,const std::vector<std::string>& atlst,Atom* src,const std::vector<Atom*>& bnd)
  {
    if(!weightedAtomSelector) throw NoAtomSelectorLoadedException();
    if(!atlst.size()) throw NoBondsAvailableException(); 
    std::vector<Atom*> chs=ff.toAtomList(atlst);
    return weightedAtomSelector->select(src,bnd,chs);
  }

}
#endif
