#ifndef INCLUDED_DNVREINFORCE
#define INCLUDED_DNVREINFORCE 1
#include "support/torchmol.h"

static int GLOBAL_TRAINMOL_COUNT=0;
namespace dnvreinforce
{
    class UnknownReinforcementAlgorithmException : public std::exception {};
    class AlreadyInitializedException : public std::exception {};

    class Policy : public mytorch::TorchModuleWrapper
    {
    public:
        Policy(const std::string& torchscriptfile)
        {
            actor=new mytorch::ConvertedTorchModule<mytorch::TorchscriptModule>(new mytorch::TorchscriptModule(torchscriptfile));
            //actor = torch::jit::load(torchscriptfile);
        }
        Policy(torch::nn::Module& tm) {actor = new mytorch::ConvertedTorchModule<torch::nn::Module>(tm);}
        Policy(mytorch::TorchModuleWrapper& tm) {actor = &tm;}

        torch::Tensor forward(mytorch::TorchTensorList x) override {return actor->forward(x[0]);}
        torch::Tensor forward(torch::Tensor x) override {return actor->forward(x);}
        std::pair<int,torch::Tensor> select_action(torch::Tensor state)
        {
            auto probs = F::softmax(this->forward(state).squeeze(1),-1);
            int c = probs.detach().multinomial(1).squeeze().item<int>();
            return make_pair(c,probs);
        }
        inline torch::nn::Module& getModelModule() {return *actor;}
        inline const torch::nn::Module& getModelModule() const {return *actor;}
        inline auto parameters(bool rec=true) {return actor->parameters(rec);}

    private:
        mytorch::TorchModuleWrapper* actor;
    };

    class ReinforcementAlgorithm
    {
    protected:
        std::vector<float> rewards;
        float gamma;
        Policy* policy;
        torch::optim::Optimizer* optimizer=nullptr;
        int mydev=0;
    protected:
        ReinforcementAlgorithm(int devid=0) {mydev=devid;}
    public:
        virtual int select_action(torch::Tensor state) =0;
        virtual inline void save_reward(float reward) {rewards.push_back(reward);}
        virtual double finish_episode(bool train=true,bool existsmol=true) = 0;
        virtual void modifyReward(float rew) {if(rewards.size()) rewards[rewards.size()-1]=rew;}
        virtual inline float* getRewardPointer() {if(rewards.size()) return &(rewards[rewards.size()-1]); else return nullptr;}
        virtual inline Policy& getPolicy() {return *policy;}
        virtual inline const Policy& getPolicy() const {return *policy;}
        virtual inline torch::optim::Optimizer* getOptimizer() {return optimizer;}

        virtual inline void replaceOptimizer(torch::optim::Optimizer* opt,bool delt=false) {if(delt) delete optimizer; optimizer = opt;}
    };

    class Reinforce : public ReinforcementAlgorithm
    {
    public:
        Reinforce(const std::string& torchscriptfile,float learning_rate=1e-3, float gamma=0.995) : Reinforce(new Policy(torchscriptfile),learning_rate,gamma) {}
        Reinforce(Policy& pol,float learning_rate=1e-3, float gamma=0.995) : Reinforce(&pol,learning_rate,gamma) {}
        
        Reinforce(Policy* pol,float learning_rate=1e-3, float gam=0.995,int devid=0) : ReinforcementAlgorithm(devid)
        {
            gamma=gam;
            policy = pol;
            optimizer = new torch::optim::Adam(policy->parameters(true), torch::optim::AdamOptions(learning_rate));
        }
        int select_action(torch::Tensor state) override;
        virtual double finish_episode(bool train=true,bool existsmol=true) override;

    private:
        std::vector<torch::Tensor> saved_log_probs;
    };

    class EnsembleReinforce
    {
      std::vector<ReinforcementAlgorithm*> reinforces;
      bool initialized=false;
    public:
        EnsembleReinforce() {}
        EnsembleReinforce(std::vector<Policy>& policies,int reinforce_type=0) : EnsembleReinforce()
        {
            this->init(policies,reinforce_type);
        }

        void init(std::vector<Policy>& policies,int reinforce_type=0)
        {
            if(initialized) throw AlreadyInitializedException();
            initialized=true;
            for(Policy& pol : policies)
            {
                if(reinforce_type==0) reinforces.push_back(new Reinforce(&pol));
                else throw UnknownReinforcementAlgorithmException();
            }
        }

        void mergeModels()
        {
            if(reinforces.size()<=1) return;
            std::vector<torch::Tensor> params;
            for(auto& parm : reinforces[0]->getOptimizer()->parameters()) params.push_back(parm.clone());
            for(int i=1;i<reinforces.size();i++)
            {
                int j=0;
                for(auto& parm : reinforces[i]->getOptimizer()->parameters()) params[j++]+=parm;
            }
            for(torch::Tensor& tens : params) tens/=(int)reinforces.size();
            for(int i=0;i<reinforces.size();i++)
            {
                int j=0;
                for(auto& parm : reinforces[i]->getOptimizer()->parameters()) parm=params[j++];
            }
            std::cout << "Compiled and merged all models\n";
        }

        inline ReinforcementAlgorithm* operator[](int idx) const {return reinforces[idx];}
    };
}

static const float FAIL_REWARD=-2.5,WIN_REWARD=1.0,ERR_REWARD=-0.05,STEP_REWARD=-0.01;

int dnvreinforce::Reinforce::select_action(torch::Tensor state) //ERR: State vector has no gradient
{
    std::pair<int,torch::Tensor> action = policy->select_action(state);
    saved_log_probs.push_back(get<1>(action)[get<0>(action)].log());
    //saved_log_probs.push_back(get<1>(action).index({get<0>(action)}).log());
    return get<0>(action);
}
static int DEBUG_COUNTER_1=0;
double dnvreinforce::Reinforce::finish_episode(bool train,bool existsmol)
{
    if(!rewards.size())
    {
        rewards.clear();
        return 0.0;
    }
    #ifdef SERIAL
    Molecule* finalmolecule=molenvmt->getMolecule();
	int osz=finalmolecule->getSize();
	int csz=molenvmt->getNodes().size();
    #else
    Molecule* finalmolecule=molenvmt[mydev]->getMolecule();
    int osz=finalmolecule->getSize(),csz=molenvmt[mydev]->getNodes().size();
    #endif
	std::cout << "Comparing mol sizes: " << csz <<" vs " << osz << "\n";

    if(train) std::cout << "At molecule "+std::to_string(GLOBAL_TRAINMOL_COUNT++)+" with existence="+to_string(existsmol)+": ";
    std::vector<torch::Tensor> policy_loss;
    float R = 0;
    for (auto r = rewards.rbegin(); r != rewards.rend(); ++r)
    {
        R = *r + gamma * R;
        policy_loss.push_back(-saved_log_probs.back() * R);
        saved_log_probs.pop_back();
        std::cout <<*r<<" ";
    }
    if(train) std::cout <<".\tFor a discounted reward of "+std::to_string(R)+"\n";
    optimizer->zero_grad();
    torch::Tensor loss = torch::stack(policy_loss).sum();
    if(train)
    {
        cout << "Backward called\n";
        loss.backward();
        #ifdef SERIAL
        for(auto& gparm : molenvmt->getMoleculeFeaturizer().parameters())
        {
            if(gparm.grad().defined()) cout <<"Gradient bounds: "<< gparm.grad().max().item<float>() << " " << gparm.grad().min().item<float>() << "\n";
            else cerr << "Error! Gradient not defined for convolver\n";
            /*{
                gparm.grad().clamp_(-1,1);
            }*/
        }
        #endif
        optimizer->step();
    }
    rewards.clear();
    return R;
}

// Language: cpp
// Path: support/reinforce.cpp
// Written with the help of GitHub Copilot
#endif
