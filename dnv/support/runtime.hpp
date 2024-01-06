#ifndef DNV_INCLUDED_RUNTIME
#define DNV_INCLUDED_RUNTIME 1
#include "graph/Molecule.hpp"

namespace runtime
{
    class InvalidInputParametersException : std::exception {};
    ForceField* runtimeff=nullptr;

    struct InputParameters
    {
        std::string fnam="",format="pdb",progname="dnv_protgen",flagcodes="0000";
        int num=25,osc=5,seednum=30,step=1,maxtries=300,ligtries=50;
        float prunefrac=0.67;
        unsigned short int gentype=0; // 0 - Unselected, 1 - Standard, 2 - Peptide
        std::vector<int> size={15},hsres;
        std::vector<std::string> seednames;
        std::string source_folder="data";
        std::string ff_file=source_folder+"/final_ff_parameters.ffin",cats_file=source_folder+"/categories.data",charges_file=source_folder+"/charges.data",defs_file=source_folder+"/definitions.data";
        std::string bondsfile=source_folder+"/itps/bondtypes.itp",anglesfile=source_folder+"/itps/angletypes.itp",dihedralsfile=source_folder+"/itps/dihedraltypes.itp",impropersfile=source_folder+"/itps/impropertypes.itp";
        std::string protein_ff=source_folder+"/protff.ffin",resid_ff=source_folder+"/protein.rtp";
        std::string extra_atypes="",extra_rules="";
        std::string extra_bonds="",extra_angles="",extra_dihs="";
        double dielectric=30,restrain_distance=0.8,spread=1.2;
        bool autospread=true,prune=false;
        std::string prunemolfile="prunesrc.pdb";

        inline void updateFFData(int gentype=1)
        {
            switch(gentype)
            {
            case 1:
                ff_file=source_folder+"/final_ff_parameters.ffin";
                break;
            case 2:
                ff_file=source_folder+"/aminoacids_dnv.ffin";
                break;
            default:
                ff_file=source_folder+"/final_ff_parameters.ffin";
            }
            cats_file=source_folder+"/categories.data";
            charges_file=source_folder+"/charges.data";
            defs_file=source_folder+"/definitions.data";
            bondsfile=source_folder+"/itps/bondtypes.itp";
            anglesfile=source_folder+"/itps/angletypes.itp";
            dihedralsfile=source_folder+"/itps/dihedraltypes.itp";
            impropersfile=source_folder+"/itps/impropertypes.itp";
            protein_ff=source_folder+"/protff.ffin";
            resid_ff=source_folder+"/protein.rtp";
        }
    };

    static std::vector<int> getResidueNumbersFromIdentifiers(const std::vector<std::string>& residIDs,int step=1)
    {
        std::vector<int> ret;
        for(std::string s : residIDs)
        {
            if(stringfx::indexOf('-',s)!=std::string::npos)
            {
            auto parts=stringfx::split(s,'-');
            int start=std::stoi(parts[0]),end=std::stoi(parts[1]);
            for(int ie=start;ie<=end;ie+=step) ret.push_back(ie);
            }
            else ret.push_back(std::stoi(s));
        }
        return ret;
    }

    static InputParameters getInputParametersFromFile(const std::string& filename) //May set global parameters
    {
        InputParameters INP;
        std::string line;
        std::ifstream loadf; loadf.open(filename);
        std::vector<std::string> residIDs;

        //All options defaults
        char s1='0',s2='0',s3='0',s4='0';
        USE_PHYSICS=false;
        QUIET=true;
        SYNTHGET=true;
        IPADDR="127.0.0.1";
        PORT=5555;
        SYNTHPENALTY=350; //Strong - Don't accept without //Default - Set very high to get exclusively synthesizable molecules (mostly)

        std::string sizes="15";

        while(getline(loadf,line))
        {
            std::vector<std::string> tokens=stringfx::split(line,'=');
            std::string idtoken=stringfx::trim(tokens[0]);
            if(idtoken=="Protein") INP.fnam=stringfx::trim(tokens[1]);
            else if(idtoken=="ProteinFormat") INP.format=stringfx::trim(tokens[1]);
            else if(idtoken=="Dielectric") INP.dielectric=std::stod(stringfx::trim(tokens[1]));
            else if(idtoken=="RestrainDistance") INP.restrain_distance=std::stod(stringfx::trim(tokens[1]));
            else if(idtoken=="SeedCount") INP.seednum=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="Step") INP.step=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="Sizes" || idtoken=="Sizes") sizes=stringfx::trim(tokens[1]);
            else if(idtoken=="MolCount") INP.num=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="SourceFolder") INP.source_folder=stringfx::trim(tokens[1]);
            else if(idtoken=="Oscillations") INP.osc=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="ProgramName") INP.progname=stringfx::trim(tokens[1]);
            else if(idtoken=="Strategy_Restrain") s1=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
            else if(idtoken=="Strategy_Superseed") s2=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
            else if(idtoken=="Strategy_Reseed") s3=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
            else if(idtoken=="Strategy_Deseed") s4=(stringfx::trim(tokens[1])=="Yes")?'1':'0';
            else if(idtoken=="Residues") residIDs=stringfx::split(tokens[1],',');
            else if(idtoken=="PruneMol") {INP.prunemolfile=stringfx::trim(tokens[1]); INP.prune=true;}
            else if(idtoken=="PruneKeepFrac") INP.prunefrac=std::stof(stringfx::trim(tokens[1]));
            else if(idtoken=="SynthCheck") SYNTHGET=(stringfx::trim(tokens[1])=="Yes");
            else if(idtoken=="SynthStrength") SYNTHPENALTY=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="Optimize") USE_PHYSICS=(stringfx::trim(tokens[1])=="Yes");
            else if(idtoken=="Quiet") QUIET=(stringfx::trim(tokens[1])=="Yes");
            else if(idtoken=="Charges") FFCHARGES=(stringfx::trim(tokens[1])=="FF");
            else if(idtoken=="IP") IPADDR=stringfx::trim(tokens[1]);
            else if(idtoken=="PORT") PORT=std::stoi(tokens[1]);
            else if(idtoken=="SeedSpread" || idtoken=="Spread") {INP.spread=std::stod(tokens[1]); INP.autospread=false;}
            else if(idtoken=="MaxTries") INP.maxtries=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="TriesPerLigand") INP.ligtries=std::stoi(stringfx::trim(tokens[1]));
            else if(idtoken=="ExtraAtoms") INP.extra_atypes+=","+stringfx::trim(tokens[1]);
            else if(idtoken=="ExtraBonds") INP.extra_bonds+=","+stringfx::trim(tokens[1]);
            else if(idtoken=="ExtraAngles") INP.extra_angles+=","+stringfx::trim(tokens[1]);
            else if(idtoken=="ExtraDihs" || idtoken=="ExtraDihedrals") INP.extra_dihs+=","+stringfx::trim(tokens[1]);
            else if(idtoken=="ExtraRules" || idtoken=="ExtraDefs" || idtoken=="ExtraDefinitions") INP.extra_rules+=","+stringfx::trim(tokens[1]);
            else if(idtoken=="Seeds")
            {
                if(stringfx::trim(tokens[1])=="All") INP.seednames={};
                else INP.seednames=stringfx::split(tokens[1],',');
            }
            else if(idtoken=="Type")
            {
                if(stringfx::trim(tokens[1])=="Peptide" || stringfx::trim(tokens[1])=="peptide" || stringfx::trim(tokens[1])=="pep" || stringfx::trim(tokens[1])=="Pep") INP.gentype=2;
                else INP.gentype=1;
            }
        }
        loadf.close();

        if(INP.gentype==0) INP.gentype=1;
        if(INP.gentype==1)
        {
            if(!INP.seednames.size()) INP.seednames={"CT1","CT2","CT3","CA","NX","CAS","NX","CPT","NY","CY","CT2x","CT1x","CP2","CE1A","CE2A","PL","SL","C3","OH1","OC"};
        }
        else if(INP.gentype==2)
        {
            if(!INP.seednames.size()) INP.seednames={"C"};
            VARIEDCOUNT=1;
            NOSEEDCHECK=1;
            ALLOW_BIPHENYL=false;
            INP.source_folder=INP.source_folder+"/aminoacids";
        }

        INP.size=getResidueNumbersFromIdentifiers(stringfx::split(sizes,','),INP.step);
        std::string strategy(1,s1); strategy+=s2; strategy+=s3; strategy+=s4;
        INP.flagcodes=strategy;
        INP.updateFFData(INP.gentype);

        INP.hsres=getResidueNumbersFromIdentifiers(residIDs);
        return INP;
    }
}

static void tuneGlobalParametersByInput(const runtime::InputParameters& INP)
{
    DIEL/=(INP.dielectric/3); //Dielectric constant is 30
    RESTRAINDISTANCE=INP.restrain_distance; //Restrain distance if enabled
}

static ForceField& loadForceFieldFromInput(const runtime::InputParameters& INP,bool verbose=true)
{
    ForceField& ff=*(new ForceField(INP.ff_file,INP.cats_file));
    ff.loadBondParameters(INP.bondsfile,INP.anglesfile,INP.dihedralsfile,INP.impropersfile);
    std::string extra_data="";
    if(INP.extra_atypes!="")
    {
        extra_data+="\tAtom types: "+INP.extra_atypes.substr(1)+"\n";
        for(const std::string& fnm : stringfx::split(INP.extra_atypes.substr(1),',')) {if(verbose) std::cout <<"Loading: "<< fnm << "\n"; ff.loadAdditionalAtomTypes(stringfx::trim(fnm));}
    }
    if(INP.extra_bonds!="")
    {
        extra_data+="\tBond Data: "+INP.extra_bonds.substr(1)+"\n";
        for(const std::string& fnm : stringfx::split(INP.extra_bonds.substr(1),',')) {if(verbose) std::cout <<"Loading: "<< fnm << "\n"; ff.loadLengths(stringfx::trim(fnm),false);}
    }
    if(INP.extra_angles!="")
    {
        extra_data+="\tAngle Data: "+INP.extra_angles.substr(1)+"\n";
        for(const std::string& fnm : stringfx::split(INP.extra_angles.substr(1),',')) {if(verbose) std::cout <<"Loading: "<< fnm << "\n"; ff.loadAngles(stringfx::trim(fnm),false);}
    }
    if(INP.extra_dihs!="")
    {
        extra_data+="\tDihedral Data: "+INP.extra_dihs.substr(1)+"\n";
        for(const std::string& fnm : stringfx::split(INP.extra_dihs.substr(1),',')) {if(verbose) std::cout <<"Loading: "<< fnm << "\n"; ff.loadDihedrals(stringfx::trim(fnm),false);}
    }
    ff.loadRules(INP.defs_file);
    if(INP.extra_rules!="")
    {
        extra_data+="\tExtra Rules: "+INP.extra_rules.substr(1)+"\n";
        for(const std::string& fnm : stringfx::split(INP.extra_rules.substr(1),',')) {if(verbose) std::cout <<"Loading: "<< fnm << "\n"; ff.loadRules(stringfx::trim(fnm),false);}
    }

    if(verbose && extra_data!="")
    {
        cout << "Extra FF data:\n";
        cout << extra_data << "\n";
    }
    runtime::runtimeff=&ff;
    return ff;
}

std::vector<Molecule*> getMoleculesFromSystems(const std::vector<System>& syses)
{
    std::vector<Molecule*> ret;
    for(const System& sys : syses) ret.push_back(sys.getDrugMolecule());
    return ret;
}
std::vector<Molecule*> getTopMolecules(std::vector<Molecule*>& mols,int take,bool deleterest=false)
{
    std::vector<Molecule*> ret;
    int K=0;
    std::vector<int> indices; std::vector<double> scores; for(Molecule* m : mols) {scores.push_back(required_energy(m)); indices.push_back(K++);}

    //Selection sort using scores
    for(int i=0;i<scores.size();i++)
    {
        int mI=i;
        double mS=scores[i];
        for(int j=i+1;j<scores.size();j++)
        {
            if(scores[j]<mS)
            {
                mS=scores[j];
                mI=j;
            }
        }
        if(mI!=i)
        {
            Molecule* tmol=mols[i];
            mols[i]=mols[mI];
            mols[mI]=tmol;

            double tsc=scores[i];
            scores[i]=scores[mI];
            scores[mI]=tsc;

            int tidx=indices[i];
            indices[i]=indices[mI];
            indices[mI]=tidx;
        }
    }
    for(int i=0;i<take;i++) ret.push_back(mols[i]);
    if(deleterest)
    {
        for(int i=take;i<mols.size();i++) delete mols[i];
    }
    return ret;
}

static void standardGeneration(const std::string& inputfile,bool verbose=true,bool prompt=true,const std::string& prefix="result_prot_", void (*global_overrides)(ForceField*) = nullptr)
{
    //Load input parameters
    runtime::InputParameters INP=runtime::getInputParametersFromFile(inputfile);

    //Tune global parameters
    tuneGlobalParametersByInput(INP);
    quickgeom::Cuboid* region=nullptr;

    //Load the forcefield
    ForceField& ff=loadForceFieldFromInput(INP,verbose);

    std::vector<int> hsres=INP.hsres;
    if(verbose)
    {
        cout << "Under: '"<<INP.progname<<"': Generating "<<INP.num<<" drugs of size="<<INP.size<<" using "<<INP.seednum<<" seeds, accumulating the best for every "<<INP.osc<<" oscillations.\n";
        cout << "Using the protein file: "<<INP.fnam<<"\tin the "<<INP.format<<" format.\n";
        if(INP.flagcodes!="0000") cout << "Running with: "<<((INP.flagcodes[0]!='0')?"restraint "+to_string(RESTRAINDISTANCE)+"A ":"")<<((INP.flagcodes[1]!='0')?"superseed ":"")<<((INP.flagcodes[2]!='0')?"reseed ":"")<<((INP.flagcodes[2]!='0')?"deseed ":"")<<"\n";
        else std::cout << "NO restraint, reseed, deseed, or superseed\n";
        std::cout << "Using dielectric: " << DIEL << "\n";
        cout << "Using residue IDs: "<<hsres<<" with seeds: "<<INP.seednames<<"\n";
    }
    if(prompt)
    {
        int tempv;
        std::cout << "Type anything to continue ... ";
        cin >> tempv;
    }

    IncompleteForceField protff(INP.protein_ff); protff.loadResidues(INP.resid_ff);
    if(INP.gentype==2) VCLIST=&(ff.getCategory("count"));

    Protein p(INP.fnam,protff,hsres,INP.format);
    cout << "Protein loaded\n";
    std::vector<Atom*> seeds;
    for(std::string s : INP.seednames) seeds.push_back(ff.getAtom(stringfx::trim(s)));
    cout << "Seeds defined\n";

    Molecule* templatemol=(INP.prune)?new Molecule(INP.prunemolfile,ff,"pdb",true):nullptr;

    if(global_overrides) global_overrides(&ff);

    for(int J : INP.size)
    {
        int K=0;
        std::cout << "Generating for size: "<<J<<"\n";
        std::vector<System> mols=p.generateLigands(ff,INP.num,J,INP.progname+"_"+to_string(J)+".log",seeds,INP.seednum,!INP.autospread,INP.spread,1,INP.osc,INP.maxtries,INP.ligtries,298,"seeds_"+INP.progname+"_"+to_string(J)+".pdb",(INP.flagcodes[0]!='0'),(INP.flagcodes[1]!='0'),(INP.flagcodes[2]!='0'),(INP.flagcodes[3]!='0'),region,templatemol,INP.prune);
        for(System& s : mols)
        {
            Molecule* m=s.getDrugMolecule();
            for(Atom* a : m->getAtoms()) cout <<a->toString()<<"	"<<ff.isSatisfied(a,m->getBondedAtoms(a))<<","<<a->isCyclized()<<"\n";
            std::string effname=(stringfx::indexOf("\%i",INP.progname)!=std::string::npos)?stringfx::replace(INP.progname,"\%i",to_string(J)):INP.progname+"_"+to_string(J);
            m->dumpMol(prefix+effname+"_"+to_string(K)+".pdb",&ff);
            m->dumpMolAsMol2(prefix+effname+"_"+to_string(K)+".mol2",&ff,false,true);
            K++;
        }
    }
}
static std::vector<std::vector<Molecule*>> geneticAlgorithm(const std::string& inputfile,int takeBest=-1,std::vector<int> sizeDrift={-1,0,1},int iterlim=100,bool verbose=true,runtime::InputParameters* INPptr=nullptr,void (*global_overrides)(ForceField*) = nullptr)
{
    //Setup return variable
    std::vector<std::vector<Molecule*>> ret;

    //Load input parameters
    runtime::InputParameters INP;
    if(!INPptr) INP=runtime::getInputParametersFromFile(inputfile);
    else
    {
        *INPptr=runtime::getInputParametersFromFile(inputfile);
        INP=*INPptr;
    }

    if(verbose) std::cout << "Taking best "<<takeBest<<" of "<<INP.num<<" generations each generation\n";

    //Check for validity of input
    if(INP.size.size()>1) {std::cerr << "Error: geneticAlgorithm expects a single size to start with. Found multiple/none.\n"; throw runtime::InvalidInputParametersException();}
    int startSize=INP.size[0];
    int effno=INP.num*sizeDrift.size();
    assert(effno>1);
    if(takeBest<=0) takeBest=::max(effno/4,1);
    assert(takeBest<effno);

    //Tune global parameters
    tuneGlobalParametersByInput(INP);
    quickgeom::Cuboid* region=nullptr;

    //Load the forcefield
    ForceField& ff=loadForceFieldFromInput(INP,verbose);

    std::vector<int> hsres=INP.hsres;
    IncompleteForceField protff(INP.protein_ff); protff.loadResidues(INP.resid_ff);
    if(INP.gentype==2) VCLIST=&(ff.getCategory("count"));

    Protein p(INP.fnam,protff,hsres,INP.format);
    cout << "Protein loaded\n";
    std::vector<Atom*> seeds;
    for(std::string s : INP.seednames) seeds.push_back(ff.getAtom(stringfx::trim(s)));
    cout << "Seeds defined\n";

    //Ignoring templatemol options on purpose

    //Global overrides
    if(global_overrides) global_overrides(&ff);

    //Starting set of molecules
    std::cout << "Generating starting set at size: "<<startSize<<"\n";
    std::vector<System> seedmols=p.generateLigands(ff,INP.num*10,startSize,INP.progname+"_"+to_string(startSize)+".log",seeds,INP.seednum,!INP.autospread,INP.spread,1,INP.osc,INP.maxtries,INP.ligtries,298,"seeds_"+INP.progname+"_"+to_string(startSize)+".pdb",(INP.flagcodes[0]!='0'),(INP.flagcodes[1]!='0'),(INP.flagcodes[2]!='0'),(INP.flagcodes[3]!='0'),region,nullptr,false);
    std::vector<Molecule*> mlist=getMoleculesFromSystems(seedmols);
    mlist=getTopMolecules(mlist,takeBest,/*deleterest=*/true);
    ret.push_back(mlist);

    std::vector<Molecule*> omols=mlist;
    for(int iterno=0;iterno<iterlim;iterno++)
    {
        std::vector<Molecule*> nmols;
        for(int J : sizeDrift)
        {
            std::vector<System> gener=p.generateLigands(ff,INP.num,J,INP.progname+"_gen"+to_string(iterno)+"_"+to_string(J)+".log",omols,2,1,INP.osc,300,50,298,(INP.flagcodes[0]!='0'),(INP.flagcodes[1]!='0'),(INP.flagcodes[2]!='0'),(INP.flagcodes[3]!='0'),region,nullptr,true,INP.prunefrac);
            std::vector<Molecule*> nmlist=getMoleculesFromSystems(gener);
            for(Molecule* m : nmlist) nmols.push_back(m);
            /*int KK=0;
            for(Molecule* savedmol : omols)
            {
                KK++;
                std::vector<System> gener=p.generateLigands(ff,INP.num,J,INP.progname+"_gen"+to_string(iterno)+"_mol"+to_string(KK)+"_"+to_string(J)+".log",seeds,INP.seednum,false,INP.spread,1,INP.osc,300,50,298,"seeds_"+INP.progname+"_"+to_string(iterno)+"_"+to_string(J)+".pdb",(INP.flagcodes[0]!='0'),(INP.flagcodes[1]!='0'),(INP.flagcodes[2]!='0'),(INP.flagcodes[3]!='0'),region,savedmol,true,INP.prunefrac);
                std::vector<Molecule*> nmlist=getMoleculesFromSystems(gener);
                for(Molecule* m : nmlist) nmols.push_back(m);
            }*/
        }
        nmols=getTopMolecules(nmols,takeBest,true);
        std::ofstream mycpt; mycpt.open("results_checkpointed_best_gen"+to_string(iterno)+".dnvcpt");
        molfx::makeCheckpoint(mycpt,nmols);
        mycpt.close();
        ret.push_back(nmols);
        omols=nmols;
    }

    delete &ff;
    return ret;
}
#endif
