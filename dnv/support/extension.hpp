#ifndef INCLUDED_EXTENSIONS
#define INCLUDED_EXTENSIONS 1
#include "graph/Molecule.hpp"

namespace extension
{
  class Extension //abstract class
  {
    std::string usage_string;
    std::string name;
    int max_args=-1,min_args=-1;

  protected:
    Extension(std::string n,std::string usage_str,int marg=-1,int Marg=-1)
    {
      name=n;
      usage_string=usage_str;
      max_args=Marg;
      min_args=marg;
    }

    virtual int runExtension(int argc,char** argv) const =0;

  public:
    virtual int run(int argc,char** argv) const
    {
      if(argc>1 && std::string(argv[1])=="-h") {this->help(); return 0;}
      else if((min_args>=0 && argc-1<min_args) || (max_args>=0 && argc-1>max_args)) {cerr << "Argument count doesn't seem to match. Check usage (below)\n"; cerr << usage_string << "\n";return 1;}
      else return this->runExtension(argc,argv);
    }

    inline virtual void help() const {cout << usage_string << "\n";}
  };
}
#endif
