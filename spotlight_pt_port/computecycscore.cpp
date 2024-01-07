#define STATICDATA 1
#define NOJSON 1
#define NOZMQ 1
#define SERIAL 1
#include "support/extension.hpp"
using extension::Extension;
//PURPOSE: To convert CHARMM FF based DeNovo molecules to SMILES Strings (From PDB)

static float fractionCyclizedScore(Molecule* mol)
{
  int nCyc=0;
  for(Atom* at : mol->getAtoms())
  {
    if(at->isCyclized()) nCyc++;
  }
  return (float)nCyc/(float)mol->getSize();
}

class ExtensionExec : public Extension
{
public:
  ExtensionExec(std::string name,std::string help,int minarg=-1,int maxarg=-1) : Extension(name,help,minarg,maxarg) {}

	int runExtension(int argc,char** argv)const override
	{
		ForceField ff(DNV_ROOT+"data/final_ff_parameters.ffin",DNV_ROOT+"/data/categories.data");
		ff.loadBondParameters();
		std::string inf=argv[1],outf=(argc>2)?argv[2]:"scores.spc";
		std::ifstream di; di.open(inf);
		std::ofstream of; of.open(outf);

		while(true)
		{
			std::string fnam; di >> fnam;
			if(di.eof()) break;
			Molecule* loadmol=new Molecule(fnam,ff,"pdb",true);
			of << fnam << " " << fractionCyclizedScore(loadmol) << "\n";
			delete loadmol;
		}
		of.close();
		cout << "Completed processing\n";
	}
};

int main(int argc,char** argv)
{
  std::string name="dnv-smigen",params="mol_list [out_file]";
  std::string extras=" \
mol_list: A file containing a list of filenames to convert. Each filename must map to a PDB file of a small molecule\n \
out_file: An output file to write the SMILES into (Default: output.smi; Also writes to stdout)";
	int min_args=1;

  std::string help_string=name+";\tUsage: "+name+" "+params+"\n"+extras;
  auto ext = ExtensionExec(name,help_string,min_args);
  return ext.run(argc,argv);
}
