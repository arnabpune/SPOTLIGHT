#ifndef INCLUDED_MYTORCH
#define INCLUDED_MYTORCH 1
#include "torch/torch.h"
#include <functional>
//#include "boost/variant.hpp"
#include <torch/script.h>
#include <memory>

using namespace torch::indexing;

using torch::nn::Module;
using torch::Tensor;

using std::placeholders::_1;
using std::placeholders::_2;

namespace mytorch
{
  typedef int64_t TorchSizeElement;
  typedef int64_t TorchIndexElement;

  std::vector<TorchSizeElement> getShape(const Tensor& t)
  {
    std::vector<TorchSizeElement> ret;
    for(int k=0;k<t.dim();k++) ret.push_back(t.size(k));
    return ret;
  }

  namespace activations
  {
    inline Tensor relu(Tensor& v) {return torch::nn::functional::relu(v,false);}
  }

  class TorchLayer : public Module
  {
    torch::nn::Module* mymodule;
    Tensor (*activation)(const Tensor&) = nullptr;
    bool activn=false;

  public:
    TorchLayer(torch::nn::Module* mymod,Tensor (*activ)(const Tensor& in)=nullptr)
    {
      mymodule=mymod;
      if(activ!=nullptr)
      {
        activn=true;
        activation=activ;
      }
    }

    inline Tensor forward(Tensor& x) 
    {
      x=mymodule->forward(x);
      if(activn) return activation(x);
      else return x;
    }

    inline auto parameters() {return mymodule->parameters();}
  };

  class TorchModel : public Module
  {
    std::vector<TorchSizeElement> inputshape,outputshape;
    //std::vector<boost::variant<torch::nn::ModuleHolder<torch::nn::LinearImpl>*,torch::nn::ModuleHolder<torch::nn::Conv2dImpl>*>> layers;
    std::vector<TorchLayer*> layers;
    //std::vector<Tensor (*)(Tensor in, ...)> acts;
    std::vector<Tensor> parameterlist;
    torch::optim::Optimizer* optim=nullptr;
  public:
    TorchModel(const std::vector<TorchSizeElement>& inshape,const std::vector<TorchSizeElement>& outshape) : Module()
    {
      inputshape=inshape;
      outputshape=outshape;
    }

    Tensor forward(Tensor& x)
    {
      for(auto* lay : this->layers)
      {
        x=lay->forward(x);
        std::cout << "\t" << getShape(x) << "\n";
      }
      return x;
    }

    torch::nn::ParameterList parameters()
    {
      torch::nn::ParameterList mylist;
      for(Tensor& t : parameterlist) mylist.get()->append(t);
      return mylist;
      //return torch::nn::ParameterList(parameterlist);
    }

    //template<class T> inline void addLayer(torch::nn::ModuleHolder<T>* mh,Tensor (*activation)(Tensor in, ...)=nullptr) {addLayer(mh->get(),activation);}

    template<typename... Args> void addLayer(TorchLayer* layer)
    {
      this->layers.push_back(layer); //This time it must be a constructed object pointer, and not just the class name
      for(auto param : this->layers[this->layers.size()-1]->parameters()) this->parameterlist.push_back(param);
      return;
    }

    std::vector<TorchSizeElement> getInputShape() {return this->inputshape;}
    std::vector<TorchSizeElement> getOutputShape() {return this->outputshape;}

    inline Tensor operator()(Tensor& x) {return this->forward(x);}
  };

  typedef std::vector<torch::Tensor> TorchTensorList;

  class TorchModuleWrapper : public torch::nn::Module
  {
  public:
    virtual torch::Tensor forward(TorchTensorList x) {return x[0];} //Overridable
    virtual torch::Tensor forward(torch::Tensor x) {return this->forward(TorchTensorList(1,x));}
    torch::Tensor operator()(torch::Tensor x) {return this->forward(x);}
  };
  template<class T> class ConvertedTorchModule : public TorchModuleWrapper
  {
    T* mymodule;
  public:
    ConvertedTorchModule(T& t) {mymodule=&t;}
    ConvertedTorchModule(T* t) {mymodule=t;}

    virtual torch::Tensor forward(TorchTensorList x) override {return mymodule->forward(x[0]);}
  };


  static bool safeTorchscriptModuleLoad(const std::string& tsfile,torch::jit::script::Module& molfeat,std::string errormsg="Error loading torchscript")
  {
    try {molfeat=torch::jit::load(tsfile); return true;}
    catch(const c10::Error& e) {std::cerr << errormsg << "\n"; return false;}
  }
  static torch::Tensor applyTSModule(torch::jit::script::Module& modl,const std::vector<torch::Tensor>& ins)
  {
    std::vector<torch::jit::IValue> inputs; for(torch::Tensor tens : ins) inputs.push_back(tens);
    return modl.forward(inputs).toTensor();
  }
  inline static torch::Tensor applyTSModule(torch::jit::script::Module& modl,const torch::Tensor& ins) {return applyTSModule(modl,std::vector<torch::Tensor>(1,ins));}

  class TorchscriptModule : public TorchModuleWrapper
  {
    torch::jit::script::Module mod;
    bool succ=false;
  public:
    TorchscriptModule(const std::string& tsfile) {succ=safeTorchscriptModuleLoad(tsfile,mod);}

    inline torch::Tensor forward(torch::Tensor x) override {return mod.forward({x}).toTensor();}
    inline torch::Tensor forward(TorchTensorList x) override
    {
      std::vector<c10::IValue> inputs; for(torch::Tensor tens : x) inputs.push_back(tens);
      return mod.forward(inputs).toTensor();
    }
    inline torch::jit::script::Module& getTorchscriptModule() {return mod;}
    inline const torch::jit::script::Module& getTorchscriptModule() const {return mod;}
  };

  namespace modules
  {
    /* A graph convolution module in PyTorch C++ */
    class GraphConvolutionImpl : public TorchModuleWrapper
    {
      torch::nn::Linear linear{nullptr};
      torch::Tensor (*activation)(const torch::Tensor&) = nullptr;
      bool activn=false;

    public:
      GraphConvolutionImpl(int64_t in_features, int64_t out_features, Tensor (*activ)(const Tensor& in)=nullptr) : TorchModuleWrapper()
      {
        linear = register_module("linear", torch::nn::Linear(in_features, out_features));
        if(activ!=nullptr)
        {
          activn=true;
          activation=activ;
        }
      }

      Tensor forward(Tensor x, Tensor adj)
      {
        Tensor output = torch::matmul(adj, linear->forward(x));
        output=output/(torch::sum(adj,1).unsqueeze(-1)+1e-10);
        if(activn) return activation(output);
        else return output;
      }

      Tensor forward(TorchTensorList x) override {return this->forward(x[0],x[1]);}
    };
    TORCH_MODULE(GraphConvolution);

    class IterativeGraphConvolutionImpl : public TorchModuleWrapper
    {
      GraphConvolution gc{nullptr};
      int64_t numiters=0;
      bool activn=false,activ_once=false;
      torch::Tensor (*activation)(const torch::Tensor&) = nullptr;

    public:
      IterativeGraphConvolutionImpl(int64_t in_features, int64_t out_features, int64_t numiters, Tensor (*activ)(const Tensor& in)=nullptr, bool activ_once=false) : TorchModuleWrapper()
      {
        gc = register_module("gc", GraphConvolution(in_features, out_features));
        this->numiters=numiters;
        if(activ!=nullptr)
        {
          activn=true;
          activation=activ;
        }
        this->activ_once=activ_once;
      }

      Tensor forward(Tensor x, Tensor adj)
      {
        Tensor output = x;
        for(int64_t i=0; i<numiters; i++)
        {
          output = gc->forward(output, adj);
          if(activn && (!activ_once || i==numiters-1)) output = activation(output);
        }
        return output;
      }

      Tensor forward(TorchTensorList x) override {return this->forward(x[0],x[1]);}
    };
    TORCH_MODULE(IterativeGraphConvolution);

    // Graph Attention Network
    class GraphAttentionImpl : public TorchModuleWrapper
    {
      //torch::nn::Linear linear{nullptr};
      std::vector<torch::nn::Linear> linears; 
      torch::Tensor (*activation)(const torch::Tensor&) = nullptr;
      bool activn=false;
      float alpha=1.0;
      float dropout=0.3;
      bool is_sparse=false,concat=false;
      int heads=1;

    public:
      GraphAttentionImpl(int64_t in_features, int64_t out_features, int heads=1, Tensor (*activ)(const Tensor& in)=nullptr, bool concat=true,float alpha=1.0, float dropout=0.0, bool is_sparse=false) : TorchModuleWrapper()
      {
        this->heads = heads;
        for(int i=0; i<heads; i++) linears.push_back(register_module("linear"+std::to_string(i), torch::nn::Linear(in_features, out_features)));
        this->alpha=alpha;
        this->dropout=dropout;
        this->is_sparse=is_sparse;
        this->concat=concat;
        if(activ!=nullptr)
        {
          activn=true;
          activation=activ;
        }
      }

      Tensor forward(Tensor x, Tensor adj)
      {
        std::vector<Tensor> outputs;
        Tensor h_prime;
        for(int i=0;i<heads;i++)
        {
          Tensor wh = linears[i]->forward(x);
          Tensor a_input = torch::matmul(wh, wh.transpose(1,2));
          Tensor e = torch::exp(a_input * -alpha);
          Tensor zero_vec = -9e15*torch::ones_like(e);
          Tensor attention = torch::where(adj > 0, e, zero_vec);
          attention = torch::softmax(attention, 1);
          attention = torch::dropout(attention, dropout, is_sparse);
          h_prime = torch::matmul(attention, wh);
          outputs.push_back(h_prime);
        }
        if(concat) h_prime = torch::cat(outputs,1);
        else h_prime = torch::mean(torch::stack(outputs),0);
        if(activn) return activation(h_prime);
        else return h_prime;
      }

      Tensor forward(TorchTensorList x) override {return this->forward(x[0],x[1]);}
    };
    TORCH_MODULE(GraphAttention);

    // Iterative Graph Attention Network
    class IterativeGraphAttentionImpl : public TorchModuleWrapper
    {
      GraphAttention ga{nullptr};
      int64_t numiters=0;
      bool activn=false,activ_once=false;
      torch::Tensor (*activation)(const torch::Tensor&) = nullptr;

    public:
      IterativeGraphAttentionImpl(int64_t in_features, int64_t out_features, int64_t numiters, int heads=1, Tensor (*activ)(const Tensor& in)=nullptr, bool activ_once=false) : TorchModuleWrapper()
      {
        ga = register_module("ga", GraphAttention(in_features, out_features, heads, nullptr, false));
        this->numiters=numiters;
        if(activ!=nullptr)
        {
          activn=true;
          activation=activ;
        }
        this->activ_once=activ_once;
      }

      Tensor forward(Tensor x, Tensor adj)
      {
        Tensor output = x;
        for(int64_t i=0; i<numiters; i++)
        {
          output = ga->forward(output, adj);
          if(activn && (!activ_once || i==numiters-1)) output = activation(output);
        }
        return output;
      }

      Tensor forward(TorchTensorList x) override {return this->forward(x[0],x[1]);}
    };
  }

  static std::vector<torch::Tensor> combineParameters(const std::vector<torch::Tensor>& a,const std::vector<torch::Tensor>& b)
  {
    std::vector<torch::Tensor> out;
    for(auto& t : a) out.push_back(t);
    for(auto& t : b) out.push_back(t);
    return out;
  }
}
#endif
