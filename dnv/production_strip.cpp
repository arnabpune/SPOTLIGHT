#define STATICDATA 1
#define NOJSON 1
#define NOINTERFACE 1
//#define USE_CUDA 1 //Disable to remove CUDA support
//#define SERIAL 1 //Enable to run serial job
//#include "graph/Molecule.hpp"
#include "support/runtime.hpp"

static inline void global_overrides(ForceField* ff=nullptr) {} //Make any manual changes to global variables here. They will only apply directly in the molecule and seed generation parts.
int main(int argc,char** argv)
{
  #ifdef USE_CUDA
  N_CURAND=(1<<24)*12;
  #endif

  std::string filename=argv[1];
  standardGeneration(filename,/*verbose=*/true,/*prompt=*/true,/*prefix=*/"result_prot_",global_overrides);
  return 0;
}
