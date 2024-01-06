#include "cudarng_mono.h"
#include <random>
#include <stdlib.h>
#include <stddef.h>
#include <iostream>

int cuRNG::MAX_RANDS = (1 << 24);
//int cuRNG::WAIT = 20000;
// Start multiple generators for parallelism
namespace cuRNG
{
    inline void resetIndices(int maxn=-1) 
    {
        for(int i=0;i<ngen;i++) index[i]=i;
        if(maxn!=-1) for(int i=0;i<ngen;i++) max_index[i]=maxn;
    }

    void init(int genc,int seed,int max_n,bool preload)
    {
        if(seed==-1)
        {
            std::random_device rd;
            seed = rd();
        }
        d_rands = (float *)malloc(max_n * sizeof(float));
        ngen=genc;
        //Initialize the generator
        curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);

        //Set the seed
        curandSetPseudoRandomGeneratorSeed(gen, seed);

        //Allocate GPU memory for the random numbers
        cudaMalloc((void **)&h_rands, max_n * sizeof(float));

        //Initialize multiple generators
        for(int i=0;i<ngen;i++)
        {
            index.push_back(i);
            max_index.push_back(0);
        }
        if(preload) generateRandomNumbers(max_n);
        initialized=true;
    }


    float* generateRandomNumbers(long n)
    {
        lock=true;
        //usleep(WAIT);
        resetIndices(n);
        curandGenerateUniform(gen, h_rands, n);
        cudaMemcpy(d_rands, h_rands, n * sizeof(float), cudaMemcpyDeviceToHost);
        lock=false;
        return d_rands;
    }

    void destroy()
    {
        while(lock) {}
        curandDestroyGenerator(gen);
        cudaFree(h_rands);
        free(d_rands);
        initialized=false;
    }

    float getRN(int devid,int refill)
    {
        //if(index[devid]%1000000==0) std::cout << std::to_string(index[devid]) + " of " + std::to_string(max_index[devid]) << "\n";
        if(!initialized) {std::cout << "cuRNG not initialized!\n"; exit(1);}
        while(lock) {}
        if(index[devid]>=max_index[devid])
        {
            if(refill==-1) refill=MAX_RANDS;
            if(!lock)
            {
                std::cout << "# Calling CUDA RNG with n="+std::to_string(refill)+"\n";
                generateRandomNumbers(refill);
            }
            while(lock) {}
        }
        float ret = d_rands[index[devid]];
        index[devid]+=ngen;
        return ret;
    }
}

/*int main(int argc,char** argv)
{
    long N = 1 << 26;
    std::cout << "N=" << N << "\n";
    cuRNG::init();
    for(int times=0;times<100;times++)
    {
        float *d_randomNumbers = cuRNG::generateRandomNumbers(N);
        std::cout << d_randomNumbers[0] << "\n";
        free(d_randomNumbers);
    }
    cuRNG::destroy();
    std::cout << "Done!\n";
}*/
