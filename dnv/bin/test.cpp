#define FFCHARGES 1
#define STATICDATA 1
#include "graph/Molecule.hpp"

int main(int argc,char** argv)
{
 DIEL/=3.0;
 quickgeom::Cuboid* region=nullptr;
 std::string fnam="",format="pdb",progname=(argc>7)?argv[7]:"dnv_protgen",flagcodes=(argc>8)?argv[8]:"0000";
 int size=(argc>3)?std::stoi(argv[3]):15,num=(argc>4)?std::stoi(argv[4]):25,osc=(argc>5)?std::stoi(argv[5]):5,seednum=(argc>6)?std::stoi(argv[6]):30,tempv;
 if(argc>1) fnam=argv[1];
 if(argc>2) format=argv[2];
 else {cout << "Enter filename: "; cin >>fnam;}
 cout << "Under: '"<<progname<<"': Generating "<<num<<" drugs of size="<<size<<" using "<<seednum<<" seeds, accumulating the best for every "<<osc<<" oscillations\n";
 if(flagcodes!="0000") cout << "Running with: "<<((flagcodes[0]!='0')?"restraint "+to_string(RESTRAINDISTANCE)+"A ":"")<<((flagcodes[1]!='0')?"superseed ":"")<<((flagcodes[2]!='0')?"reseed ":"")<<((flagcodes[2]!='0')?"deseed ":"")<<"\n";
 cin >> tempv;
 ForceField ff("/home/venkata/dnv/data/final_ff_parameters.ffin","/home/venkata/dnv/data/categories.data");
 ff.loadBondParameters();
 ff.loadRules("/home/venkata/dnvs/data/definitions.data");
 IncompleteForceField protff("/home/venkata/dnv/data/protff.ffin"); protff.loadResidues("/home/venkata/dnv/data/protein.rtp");
sad
 std::vector<int> hsres={800,900};

 Protein p(fnam,protff,hsres,format);
 std::vector<Atom*> seeds;
 seeds.push_back(ff.getAtom("CT1"));
 seeds.push_back(ff.getAtom("CT2"));
 seeds.push_back(ff.getAtom("CT3"));
 seeds.push_back(ff.getAtom("CA"));
 seeds.push_back(ff.getAtom("NX"));
 seeds.push_back(ff.getAtom("CAS"));
 seeds.push_back(ff.getAtom("NX"));
 seeds.push_back(ff.getAtom("CPT"));
 seeds.push_back(ff.getAtom("NY"));
 seeds.push_back(ff.getAtom("CY"));
 seeds.push_back(ff.getAtom("CT2x"));
 seeds.push_back(ff.getAtom("CT1x"));
 seeds.push_back(ff.getAtom("CP2"));
 seeds.push_back(ff.getAtom("CE1A"));
 seeds.push_back(ff.getAtom("CE2A"));
 seeds.push_back(ff.getAtom("PL"));
 seeds.push_back(ff.getAtom("SL"));
 seeds.push_back(ff.getAtom("C3"));
 seeds.push_back(ff.getAtom("OH1"));
 seeds.push_back(ff.getAtom("OC"));


 Molecule* templatemol=(argc>9)?new Molecule(std::string(argv[9]),ff,"pdb",true):nullptr;
 int K=0;
 std::vector<System> mols=p.generateLigands(ff,num,size,progname+".log",seeds,seednum,false,1.2,1,osc,300,50,298,"seeds_"+progname+".pdb",(flagcodes[0]!='0'),(flagcodes[1]!='0'),(flagcodes[2]!='0'),(flagcodes[3]!='0'),region,templatemol,false);
 for(System& s : mols)
 {
  Molecule* m=s.getDrugMolecule();
  for(Atom* a : m->getAtoms()) cout <<a->toString()<<"	"<<ff.isSatisfied(a,m->getBondedAtoms(a))<<","<<a->isCyclized()<<"\n";
  m->dumpMol("result_prot_"+progname+"_"+to_string(K)+".pdb",&ff);
  K++;
 }
}

