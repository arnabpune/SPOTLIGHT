#pragma once
#include <cuda.h>
#include <curand.h>
#include <vector>

#include "driver_types.h"
#include <cuda_runtime.h>

namespace cuRNG
{
    extern int MAX_RANDS;
    static std::vector<curandGenerator_t> generators;
    static float** d_rands;
    static float** gpu_randomNumbers;
    static std::vector<int> index,max_index;
    static bool initialized = false;
    void init(int ngen=1,int seed=-1,int max_n=MAX_RANDS,bool preload=true);
    float* generateRandomNumbers(long n,int devid=0);
    void destroy();
    float getRN(int devid=0,int refill=-1);
}