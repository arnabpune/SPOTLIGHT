#ifndef HPP_DNVMACROS
#define HPP_DNVMACROS 1
//#define USE_PHYSICS 1
//#define QUIET 1
#define MACROMODE 1
#include "graph/Molecule.hpp"
#define FUNCCHAR '@'
#define FUNCCHAR2 '~' //Reference functions (that use only variable name, not value)
#define VARCHAR '$'
#define MODCHAR '%'
namespace macros
{
  class BadTypeException : public exception {};
  class UndefinedEnvironmentParameterException : public exception {};
  class NoSuchVariableException : public exception {};
  class BadCallToUncallableVariableException : public exception {};
  class MacroSyntaxError : public exception {};
  class InvalidReferenceToUnsubscriptableVariableException : public exception {};
  class NoSuchParameterToObjectException : public exception {};
  class BadSubsriptIndexException : public exception {};
  class BadIterationObjectException : public exception {};
  class BadFunctionInputException : public exception {};
  class NoSuchFunctionException : public exception {};
  class InvalidMathOperationException : public exception {};
  class RetrievingExistingVariableException : public exception {};
  class NoForceFieldLoadedException : public UndefinedEnvironmentParameterException {};

  template<class T> static inline double quicknorm(const std::vector<T>& v)
  {
    double r=0;
    for(const T& i : v) r+=i*i;
    return sqrt(r);
  }
  class Variable
  {
    void* data=nullptr;
    int type;
    std::string name,typestr="";
    const ForceField* myff=nullptr;
    bool callable=false,subscriptable=false,iterable=false;
    std::vector<std::string> myStructTypes=std::vector<std::string>();

    static std::vector<std::string> typenames;
    static std::vector<std::string> structlist;
    static std::vector<std::vector<std::string>> structdef;

  public:
    Variable() {}
    Variable(const std::string& n,int tID,const std::string& dtext="",const ForceField* nff=nullptr)
    {
      name=n;
      myff=nff;
      type=tID;
      assign(dtext);
    }
    Variable(const std::string& n,int tID,void* ready,const ForceField* nff=nullptr)
    {
      name=n;
      myff=nff;
      type=tID;
      data=ready;
    }
    Variable(const std::string& n,const std::string& tstr,const std::string& dtext="",const ForceField* nff=nullptr)
    {
        name=n;
        myff=nff;
        setType(tstr);
        assign(dtext);
    }
    Variable(const std::string& n,const std::string& tstr,void* ready,const ForceField* nff=nullptr)
    {
        name=n;
        myff=nff;
        setType(tstr);
        data=ready;
    }
    /*Variable(const Variable& var)
    {

    }
    ~Variable() {clearData();}*/

    inline static void addStruct(const std::string& sn,const std::string& sref)
    {
      structlist.push_back(stringfx::trim(sn));
      std::vector<std::string> strs=stringfx::split(sref,',');
      for(std::string& s : strs) s=stringfx::trim(s);
      structdef.push_back(strs);
    }
    inline static std::vector<std::string> getStructDef(const std::string& sn)
    {
      for(int i=0;i<structlist.size();i++)
      {
        if(structlist[i]==sn) return structdef[i];
      }
      return std::vector<std::string>();
    }
    inline const std::string& getName() const {return name;}
    inline bool isCallable() const {return callable;}
    inline bool isSubscriptable() const {return subscriptable;}
    inline bool isIterable() const {return iterable;}
    inline const void* getRawObjectPointer() const {return data;}
    inline void* getRawObjectPointer() {return data;}
    inline int getType() const {return type;}
    inline int getTypeID() const {return type;}
    inline const ForceField* getForceField() const {return myff;}
    inline const std::string& getTypeString() const {return typestr;}
    inline std::string getReplacementString() const;
    void clearData();

    static int getTypeUID(const std::string& ts)
    {
      for(int i=0;i<typenames.size()-1;i++)
      {
        if(typenames[i]==ts) return i;
      }
      if(contains(structlist,ts)) return 10;
      else return -1;
    }

    void setType(const std::string& ts)
    {
      if(ts=="int") type=0;
      else if(ts=="double" || ts=="float") type=1;
      else if(ts=="string") {type=2; subscriptable=true; callable=true;}
      else if(ts=="atom") {type=3; callable=true;}
      else if(ts=="molecule") {type=4; callable=true; subscriptable=true; iterable=true;}
      else if(ts=="intvector") {type=5; subscriptable=true; callable=true; iterable=true;}
      else if(ts=="doublevector") {type=6; subscriptable=true; callable=true; iterable=true;}
      else if(ts=="stringvector") {type=7; subscriptable=true; callable=true;iterable=true;}
      else if(ts=="molvector") {type=8; subscriptable=true; callable=true; iterable=true;}
      else
      {
        if(ts.find("array")==0)
        {
          typestr=ts.substr(6,ts.find_last_of('>')-6);
          type=9;
          subscriptable=true; iterable=true; callable=true;
        }
        else
        {
          std::vector<string> nsdef=getStructDef(ts);
          if(!nsdef.size())
          {
            cerr << "Unknown type: "<<ts<<"\n";
            throw BadTypeException();
          }
          else
          {
            type=10; subscriptable=true;
            for(const std::string& st : nsdef) myStructTypes.push_back(st);
          }
        }
      }
    }
    void assign(const std::string& ls);


    inline std::string getTypedReplacementString() const  {return typenames[type]+"~"+getReplacementString();}
  };
  std::vector<std::string> Variable::typenames={"int","double","string","atom","molecule","intvector","doublevector","stringvector","molvector","array","struct"};
  std::vector<std::string> Variable::structlist={};
  std::vector<std::vector<std::string>> Variable::structdef={};

  namespace callables
  {
    inline static Variable callVar(const Variable& v,const std::string& obj,const ForceField* subff=nullptr)
    {
      double* rd=nullptr;
      switch(v.getType())
      {
        case 2:
          if(obj=="length") {return Variable(obj,"int",to_string(((std::string*)v.getRawObjectPointer())->length()),subff);}
        //Molecule
        case 3:
          if(obj=="resnum") {return Variable(obj,"int",to_string(((Atom*)v.getRawObjectPointer())->getResidueNumber()),subff);}
          else if(obj=="name") {return Variable(obj,"string",((Atom*)v.getRawObjectPointer())->toString(),subff);}
          else if(obj=="weight") {return Variable(obj,"double",to_string(((Atom*)v.getRawObjectPointer())->seek_mass()),subff);}
        case 4:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((Molecule*)v.getRawObjectPointer())->getSize()),subff);}
          else if(obj=="effectivelength" || obj=="effectivesize") {return Variable(obj,"int",to_string(((Molecule*)v.getRawObjectPointer())->getEffectiveSize()),subff);}
          else if(obj=="COG" || obj=="CentreOfGeometry")
          {
            //Eigen::Vector3d* myvec=new Eigen::Vector3d(((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry());
            std::vector<double>* retv=new std::vector<double>(3); for(int it=0;it<3;it++) retv->operator[](it)=((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry(false)(it);
            return Variable(obj,"doublevector",retv,subff);
          }
          else if(obj=="COGH" || obj=="CentreOfGeometryReduced")
          {
            //Eigen::Vector3d* myvec=new Eigen::Vector3d(((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry());
            std::vector<double>* retv=new std::vector<double>(3); for(int it=0;it<3;it++) retv->operator[](it)=((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry(true)(it);
            return Variable(obj,"doublevector",retv,subff);
          }
          else if(obj=="COM" || obj=="CentreOfMass")
          {
            //Eigen::Vector3d* myvec=new Eigen::Vector3d(((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry());
            std::vector<double>* retv=new std::vector<double>(3); for(int it=0;it<3;it++) retv->operator[](it)=((Molecule*)v.getRawObjectPointer())->getStoredCentreOfMass(false)(it);
            return Variable(obj,"doublevector",retv,subff);
          }
          else if(obj=="COMH" || obj=="CentreOfMassReduced")
          {
            //Eigen::Vector3d* myvec=new Eigen::Vector3d(((Molecule*)v.getRawObjectPointer())->getStoredCentreOfGeometry());
            std::vector<double>* retv=new std::vector<double>(3); for(int it=0;it<3;it++) retv->operator[](it)=((Molecule*)v.getRawObjectPointer())->getStoredCentreOfMass(true)(it);
            return Variable(obj,"doublevector",retv,subff);
          }
        case 5:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((std::vector<int>*)v.getRawObjectPointer())->size()),subff);}
          else if(obj=="norm" || obj=="magnitude")
          {
            rd=new double(quicknorm(*(std::vector<int>*)v.getRawObjectPointer()));
            return Variable(obj,"double",rd,subff);
          }
        case 6:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((std::vector<double>*)v.getRawObjectPointer())->size()),subff);}
          else if(obj=="norm" || obj=="magnitude")
          {
            rd=new double(quicknorm(*(std::vector<double>*)v.getRawObjectPointer()));
            return Variable(obj,"double",rd,subff);
          }
        case 7:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((std::vector<std::string>*)v.getRawObjectPointer())->size()),subff);}
        case 8:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((std::vector<Molecule*>*)v.getRawObjectPointer())->size()),subff);}
        case 9:
          if(obj=="length" || obj=="size") {return Variable(obj,"int",to_string(((std::vector<Variable>*)v.getRawObjectPointer())->size()),subff);}
      }
      throw NoSuchParameterToObjectException();
    }
    inline static Variable subscriptVar(const Variable& v,const std::string& obj,const ForceField* subff=nullptr)
    {
      switch(v.getType())
      {
        //Molecule
        case 2:
          return Variable("string_"+obj,"string",new std::string(1,((std::string*)v.getRawObjectPointer())->operator[](std::stoi(obj))),subff);
        case 4:
          try {return Variable("atom_"+obj,"atom",(Atom*)((Molecule*)v.getRawObjectPointer())->getAtoms()[std::stoi(obj)],subff);}
          catch(std::exception e) {return Variable("atom_"+obj,"atom",(Atom*)((Molecule*)v.getRawObjectPointer())->getAtoms()[((Molecule*)v.getRawObjectPointer())->indexOf(obj)],subff);}
        case 5:
          return Variable("int_"+obj,"int",(int*)&((std::vector<int>*)v.getRawObjectPointer())->operator[](std::stoi(obj)),subff);
        case 6:
          return Variable("int_"+obj,"int",(int*)&((std::vector<double>*)v.getRawObjectPointer())->operator[](std::stoi(obj)),subff);
        case 7:
          return Variable("int_"+obj,"int",(int*)&((std::vector<std::string>*)v.getRawObjectPointer())->operator[](std::stoi(obj)),subff);
        case 8:
          return Variable("mol_"+obj,"molecule",(Molecule*)((std::vector<Molecule*>*)v.getRawObjectPointer())->operator[](std::stoi(obj)),subff);
        case 9:
        case 10:
          return (Variable)((std::vector<Variable>*)v.getRawObjectPointer())->operator[](std::stoi(obj));
      }
      throw BadSubsriptIndexException();
    }
  }

  class Macro;
  class MacroFunction
  {
    std::string name;
  public:
    MacroFunction(const std::string& n) {name=n;}
    virtual Variable evaluate(std::string in,const Macro& mac)=0;
    inline const std::string& getName() const {return name;}
  };
  class Macro
  {
    std::vector<MacroFunction*> funcs;
    std::vector<Variable> vars;
    std::vector<std::pair<std::string,std::string>> lines;
    std::vector<int> loopstarts,loopends,whilestarts;
    std::vector<std::string> loopvars;
    std::vector<std::vector<std::string>> loopcontexts; std::vector<int> loopinds;
    std::vector<bool> ifsmatched;
    std::vector<std::pair<std::string,int>> labels;
    ForceField* myff=nullptr;
    bool procS=true;
    int layer=0;
    const Macro* linkedmac=nullptr;
  public:
    Macro() {}
    Macro(std::string macfile,bool silent=false) {loadMacro(macfile,silent);}
    Macro(const std::string& macfile, const Macro& src) : Macro(macfile,true) {linkedmac=&src;}

    void loadMacro(const std::string& file,bool silent=false)
    {
      std::ifstream ifs; ifs.open(file);
      if(!ifs.is_open()) {cerr << "File not found: "<<file<<"\n"; throw FileNotFoundException();}
      lines=std::vector<std::pair<std::string,std::string>>();
      while(!ifs.eof())
      {
        std::string cl; getline(ifs,cl);
        if(!silent) cout << cl << "\n";
        cl=stringfx::trim(cl);
        if(cl.length())
        {
          std::string comm,para;
          std::stringstream ss(cl);
          ss >> comm; getline(ss,para);
          lines.push_back(make_pair(comm,para));
        }
      }
      if(!silent) cout << "Loaded Macro: "<< lines.size() << " lines in macro\n";
      ifs.close();
    }

    Variable retreiveVar(const std::string& varname) const
    {
      if(hasVar(varname)) return getVar(varname);
      else
      {
        if(linkedmac) {return linkedmac->retreiveVar(varname);}
        else {cerr << "Retrieve failed the entire stack-search. Couldn't find '"<<varname<<"' in entire traceback \n" ; throw NoSuchVariableException();}
      }
    }
    void varLine(std::string vn,const std::string& asn)
    {
      vn= vn.substr(0,vn.length()-1);
      std::stringstream ss(asn); std::string tp, dat;
      if(hasVar(vn)) {getVar(vn).assign(asn); return;}
      ss >> tp;
      dat=stringfx::drain(ss);
      //cout << tp << "\t"<<dat<<"\n";
      vars.push_back(Variable(vn,tp,dat,myff));
    }
    std::string checkFile(const std::string& fn)
    {
      if(!fileavailable(fn))
      {
        std::string dat=ff_folder+"/../"+fn;
        if(!fileavailable(dat)) {cerr <<"Tried '"<<fn<<"' and '"<<dat<<"'. File does not exist or can't be opened\n"; throw FileNotFoundException();}
        else return dat;
      }
      else return fn;
    }
    void useData(const std::string& s)
    {
      std::stringstream ss(s); std::string tp, dat;
      ss >> tp;
      dat=stringfx::trim(stringfx::drain(ss));
      dat=checkFile(dat);
      if(tp=="forcefield")
      {
        if(myff) delete myff;
        myff=new ForceField(dat);
      }
      else if(tp=="incompleteforcefield")
      {
        if(myff) delete myff;
        myff=new IncompleteForceField(dat);
      }
      else if(tp=="categories")
      {
        if(myff) myff->loadCategories(dat);
        else {cerr << "ERR: Loading categories before forcefield.\n"; throw UndefinedEnvironmentParameterException();}
      }
      else
      {
        cerr << "Unkown file-type '"<<tp<<"' for \"using\" statement. File-type must be in ([incomplete]forcefield,categories)\n";
        throw BadFunctionInputException();
      }
    }
    void loadFFData(const std::string& inl)
    {
      if(!myff) {cerr << "ERR: Requested loading FF Data without choosing core forcefield. Did you forget the 'use' command?\n"; throw NoForceFieldLoadedException();}
      std::vector<std::string> opts=stringfx::split(inl,' ');
      if(opts[0]=="residues")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadResidues(opts[1]);
      }
      //Note: The existing data is never overwritten when loading bond/angle/dihedral parameter data. If needed, reload entire forcefield.
      if(opts[0]=="bonds")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadLengths(opts[1],false);
      }
      else if(opts[0]=="angles")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadAngles(opts[1],false);
      }
      else if(opts[0]=="dihedrals")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadDihedrals(opts[1],false);
      }
      else if(opts[0]=="impropers")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadImpropers(opts[1]);
      }
      else if(opts[0]=="rules" || opts[0]=="definitions")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadRules(opts[1]);
      }
      else if(opts[0]=="moreatoms" || opts[0]=="atomtypes")
      {
        opts[1]=checkFile(opts[1]);
        myff->loadAdditionalAtomTypes(opts[1]);
      }
    }
    inline void addFunction(MacroFunction* fx) {funcs.push_back(fx);}
    bool removeVar(const std::string& vn)
    {
      for(int i=0;i<vars.size();i++)
      {
        if(vars[i].getName()==vn)
        {
          vars.erase(vars.begin()+i);
          return true;
        }
      }
      cerr << "WARN: Cannot erase variable from varlist: varname='"<<vn<<"'\n";
      return false;
    }
    MacroFunction* getFunction(const std::string& s)
    {
      for(MacroFunction* p : funcs) {if(p->getName() == s) return p;}
      return nullptr;
    }
    inline bool hasVar(const std::string& vn) const
    {
      for(const Variable& v : vars) {if(v.getName()==vn) return true;}
      return false;
    }
    inline const Variable& getVar(const std::string& vn) const
    {
      for(const Variable& v : vars) {if(v.getName()==vn) return v;}
      cerr << "Requested variable: "<<vn<<" is not available\n";
      throw NoSuchVariableException();
    }
    inline Variable& getVar(const std::string& vn)
    {
      for(Variable& v : vars) {if(v.getName()==vn) return v;}
      cerr << "Requested variable: "<<vn<<" is not available\n";
      throw NoSuchVariableException();
    }
    inline std::string keyreplace(std::string s) const
    {
      while(s.find("NEWLINE")!=std::string::npos) s=stringfx::replace(s,"NEWLINE","\n");
      while(s.find("TAB")!=std::string::npos) s=stringfx::replace(s,"TAB","\t");
      while(s.find("SPACE")!=std::string::npos) s=stringfx::replace(s,"SPACE"," ");
      return s;
    }
    inline const ForceField* getMacroFF() const {return myff;}
    std::string processStatement(const std::string& ins,bool commside=false)
    {
      std::string ds=ins;
      int di=stringfx::lastIndexOf(VARCHAR,ds);
      if(di==-1) di=stringfx::lastIndexOf(FUNCCHAR,ds);
      if(di==-1) di=stringfx::lastIndexOf(FUNCCHAR2,ds);
      char tc;
      std::string repl="",offset="",templs(1,ds[di]);
      while(di!=-1)
      {
        templs=std::string(1,ds[di]);
        bool dol=(ds[di]==VARCHAR);
        std::string reps="",repl="";
        tc=ds[++di];
        while((tc>='A' && tc <='Z') || (tc>='a' && tc<='z') || (tc>='0' && tc<='9')) {reps+=tc; templs+=tc; tc=ds[++di];}
        if(tc=='(' && !dol) //Function
        {
          MacroFunction* fx=getFunction(reps);
          if(!fx) {cerr << "Function '"<<reps<<"' is undefined in this context\n"; throw NoSuchFunctionException();}
          templs+=tc;
          offset="";
          int brit=1;
          while(tc!=')' || brit)
          {
            tc=ds[++di];
            templs+=tc;
            if(tc==')') {brit--; if(!brit) break;}
            if(tc=='(') brit++;
            offset+=tc;
          }
          Variable nrepl=fx->evaluate(offset,*this);
          repl=nrepl.getReplacementString();
          nrepl.clearData();
        }
        else //Variable
        {
          Variable var=getVar(reps);
          while(tc=='.' || tc=='[')
          {
            templs+=tc;
            bool dot=(tc=='.');
            offset="";
            tc=ds[++di];
            while((tc>='A' && tc <='Z') || (tc>='a' && tc<='z') || (tc>='0' && tc<='9')) {offset+=tc; templs+=tc; tc=ds[++di];}
            if(!offset.length()) {cerr << ins <<":\t"<< "Trailing symbols from {.,[,]} in statement. Is it correct?\n"; throw MacroSyntaxError();}
            if(tc==']') {tc=ds[++di]; templs+=']';}
            if(dot)
            {
              if(!var.isCallable())
              {
                cerr << "Invalid 'call' to 'uncallable' variable in: '"<< ins << "'"<<" at variable: "<<var.getTypedReplacementString()<<"\n";
                throw BadCallToUncallableVariableException();
              }
              var=callables::callVar(var,offset,var.getForceField());
            }
            else
            {
              if(!var.isSubscriptable())
              {
                cerr << "Invalid 'subscript' to 'unsubscriptable' variable in: '"<< ins << "'"<<" at variable: "<<var.getTypedReplacementString()<<"\n";
                throw InvalidReferenceToUnsubscriptableVariableException();
              }
              var=callables::subscriptVar(var,offset,var.getForceField());
            }
          }
          repl=var.getReplacementString();
        }
        //cout << templs << "," << repl << "\n";
        ds=stringfx::replace(ds,templs,repl);
        di=stringfx::lastIndexOf(VARCHAR,ds);
        if(di==-1) di=stringfx::lastIndexOf(FUNCCHAR,ds);
        if(di==-1) di=stringfx::lastIndexOf(FUNCCHAR2,ds);
      }
      return keyreplace(ds);
    }
    void addStructure(const std::string& str) const
    {
      std::string n,r;
      std::stringstream ss(str); ss >> n >> r;
      cout << "Adding structure: "<<n<<" with ref: "<<r << "\n";
      Variable::addStruct(n,r);
    }
    void beginIteration(const std::string& nxt,int lno)
    {
      std::string itp,itv,pass,itt;
      std::stringstream ss(nxt); ss >> itp>>itv >> pass;
      itt=stringfx::drain(ss);
      int bsi=stringfx::indexOf('(',itt);
      if(bsi==-1) throw BadIterationObjectException();
      itt=itt.substr(bsi+1);
      int bei=stringfx::lastIndexOf(')',itt);
      if(bei==-1) throw BadIterationObjectException();
      itt=itt.substr(0,bei);
      std::vector<std::string> itercontext=stringfx::split(itt,',',false,'(',')');
      //cout << itercontext << "\n";
      loopstarts.push_back(lno);
      loopcontexts.push_back(itercontext);
      //cout << itercontext.size() << "\n";
      loopinds.push_back(0);
      loopvars.push_back(itv);
      varLine(itv+":",itp+" "+itercontext[0]);
    }
    inline void endIteration(int& lno,const std::string& end="iteration")
    {
      int reqr=loopcontexts.size()-1;
      loopinds[reqr]++;
      if(loopcontexts[reqr].size()<=loopinds[reqr])
      {
        loopcontexts.erase(loopcontexts.begin()+reqr);
        loopinds.erase(loopinds.begin()+reqr);
        loopstarts.erase(loopstarts.begin()+reqr);
        std::string vnam=loopvars[reqr];
        this->removeVar(vnam);
        loopvars.erase(loopvars.begin()+reqr);
      }
      else
      {
        getVar(loopvars[reqr]).assign(loopcontexts[reqr][loopinds[reqr]]);
        lno=loopstarts[reqr];
      }
    }
    void delVars(const std::string& ins)
    {
      std::vector<std::string> vlist=stringfx::split(stringfx::trim(ins),',');
      for(std::string vn : vlist) {if(!removeVar(vn)) throw NoSuchVariableException();}
    }
    void enterIfWhile(std::string comm,std::string cond,bool& procS,int& layer,int& lno)
    {
      if(comm=="if")
      {
        if(!procS) {layer++; return;}
        std::string exe=std::string(1,FUNCCHAR)+"matheval("+stringfx::trim(cond)+")";
        std::string ref=processStatement(exe,false);
        bool v=(int)std::stod(ref);
        if(v) {ifsmatched.push_back(true); procS=true;}
        else {ifsmatched.push_back(false); procS=false;}
      }
      else if(comm=="elif")
      {
        if(layer) return;
        bool mch=ifsmatched[ifsmatched.size()-1];
        if(mch) {procS=false; return;}
        std::string exe=std::string(1,FUNCCHAR)+"matheval("+stringfx::trim(cond)+")";
        std::string ref=processStatement(exe,false);
        bool v=(int)std::stod(ref);
        if(v) {ifsmatched[ifsmatched.size()-1]=true; procS=true;}
        else procS=false;
      }
      else if(comm=="else")
      {
        if(layer) return;
        bool mch=ifsmatched[ifsmatched.size()-1];
        if(mch) procS=false;
        else {ifsmatched[ifsmatched.size()-1]=true; procS=true;}
        return;
      }
      else if(comm=="end")
      {
        if(cond=="if")
        {
          if(layer) layer--;
          else ifsmatched.erase(ifsmatched.begin()+ifsmatched.size()-1);
        }
        else //end while
        {
          if(layer) layer--;
          else
          {
            if(procS)
            {
              int mvl=whilestarts[whilestarts.size()-1];
              whilestarts.erase(whilestarts.begin()+whilestarts.size()-1);
              lno=mvl-1;
            }
          }
        }
        if(!layer) procS=true;
      }
      else if(comm=="while")
      {
        if(!procS) {layer++; return;}
        std::string exe=std::string(1,FUNCCHAR)+"matheval("+stringfx::trim(cond)+")";
        std::string ref=processStatement(exe,false);
        bool v=(int)std::stod(ref);
        if(v) {procS=true; whilestarts.push_back(lno);}
        else procS=false;
      }
    }
    void receiveVar(const std::string& inl)
    {
      std::vector<std::string> vss=stringfx::split(inl,',',true);
      for(const std::string& st : vss)
      {
        if(hasVar(st)) {cerr << "Variable name: '"<<st<<"' (requested for retrieval) is already taken in target macro. \n"; throw RetrievingExistingVariableException();}
        Variable temp=retreiveVar(st);
        vars.push_back(temp);
      }
    }
    void callMacro(const std::string& cm);
    bool addLine(const std::string& ln)
    {
      std::string cl=stringfx::trim(ln);
      if(cl.length())
      {
        std::string comm,para;
        std::stringstream ss(cl);
        ss >> comm; getline(ss,para);
        lines.push_back(make_pair(comm,para));
        return true;
      }
      return false;
    }
    inline bool addLineandRun(const std::string& ln)
    {
      if(!addLine(ln)) return true;
      int l=lines.size()-1;
      bool succ;
      while(l<lines.size())
      {
        std::pair<std::string,std::string> rp=lines[l];
        succ=execLine(get<0>(rp),get<1>(rp),l);
        l++;
      }
      return succ;
    }
    bool execLine(const std::string& il1,const std::string& il2,int& l)
    {
      std::string comm,ref;
      comm=stringfx::trim(il1);
      if(comm[0]=='#') return true;
      ref=stringfx::trim(il2);
      comm=processStatement(comm,true);
      if(comm=="if" || comm=="elif" || comm=="else" || comm=="end" || comm=="while") enterIfWhile(comm,ref,procS,layer,l);
      if(!procS) return true;
      ref=processStatement(ref,false);
      if(comm=="using") useData(ref);
      else if(comm=="print") cout << ref;
      else if(comm=="iterate") beginIteration(ref,l);
      else if(comm=="next") endIteration(l,ref);
      else if(stringfx::indexOf("ff.",comm)==0) loadFFData(comm.substr(3)+" "+ref);
      else if(comm=="delete" || comm=="del") delVars(ref);
      else if(comm=="struct" || comm=="structure") addStructure(ref);
      else if(comm=="receive") receiveVar(ref);
      else if(comm=="call") callMacro(ref);
      else if(comm=="exit") return false;
      if(stringfx::indexOf(':',comm)!=-1) varLine(comm,ref);
      return true;
    }
    inline void execMacro() {for(int l=0;l<lines.size();l++) {if(!execLine(get<0>(lines[l]),get<1>(lines[l]),l)) break;}}
  };
  namespace inbuiltfuncs
  {
    static inline std::string traceString(const std::string& v,int bI,bool back=false,bool floats=true)
    {
      std::string ret="";
      bool st=false;
      while(bI>=0)
      {
        if(v[bI]=='\"') {if(!st) st=true; else break;}
        if(back) ret=v[bI--]+ret;
        else ret+=v[bI++];
      }
      if(back) return "\""+ret;
      else return ret+"\"";
    }
    static inline std::string traceNumber(const std::string& v,int bI,bool back=false,bool floats=true)
    {
      std::string ret="";
      while(bI>=0 && ((v[bI]>='0' && v[bI]<='9') || v[bI]=='-' || v[bI]==' ' || v[bI]=='\t' || (floats && v[bI]=='.')))
      {
        if(back) ret=v[bI--]+ret;
        else ret+=v[bI++];
      }
      return ret;
    }
    class ReferenceMacroFunction : public MacroFunction
    {
    public:
      ReferenceMacroFunction(const std::string& n) : MacroFunction(n) {}

      std::pair<std::vector<double>,std::vector<double>> readTwoVectors(const std::string& in,const Macro& mc,std::string oper="",int rs=0)
      {
        std::vector<double> r1,r2;
        std::vector<std::string> vecs=stringfx::split(in,',',false,'(',')');
        if(mc.hasVar(vecs[0]))
        {
          Variable vr=mc.getVar(vecs[0]);
          if(vr.getType()==6) r1=*((std::vector<double>*)vr.getRawObjectPointer());
          else
          {
            std::vector<int>* nv=((std::vector<int>*)vr.getRawObjectPointer());
            for(int& n : *nv) {r1.push_back(n);}
          }
        }
        else
        {
          std::vector<std::string> nums=stringfx::split(vecs[0].substr(1),',');
          for(const std::string& s : nums) r1.push_back(std::stod(s));
        }
        if(mc.hasVar(vecs[1]))
        {
          Variable vr=mc.getVar(vecs[1]);
          if(vr.getType()==6) r2=*((std::vector<double>*)vr.getRawObjectPointer());
          else
          {
            std::vector<int>* nv=((std::vector<int>*)vr.getRawObjectPointer());
            for(int& n : *nv) {r2.push_back(n);}
          }
        }
        else
        {
          std::vector<std::string> nums=stringfx::split(vecs[1].substr(1),',');
          for(const std::string& s : nums) r2.push_back(std::stod(s));
        }
        if(r1.size()!=r2.size() || (rs && r1.size()!=rs)) {cerr << "Vectors received: "<<vecs[0]<<"\t"<<vecs[1]<<": do not have same dimension for "<<oper<<"\n"; throw InvalidMathOperationException();}
        return make_pair(r1,r2);
      }
      std::vector<std::vector<double>> readManyVectors(const std::string& in,const Macro& mc,std::string oper="",int rs=0)
      {
        std::vector<std::vector<double>> ret;
        std::vector<std::string> vecs=stringfx::split(in,',',false,'(',')');
        for(int i=0;i<vecs.size();i++)
        {
          if(mc.hasVar(vecs[i]))
          {
            Variable vr=mc.getVar(vecs[i]);
            if(vr.getType()==6) ret.push_back(*((std::vector<double>*)vr.getRawObjectPointer()));
            else
            {
              std::vector<int>* nv=((std::vector<int>*)vr.getRawObjectPointer());
              std::vector<double> r1;
              for(int& n : *nv) {r1.push_back(n);}
              ret.push_back(r1);
            }
          }
        }
        return ret;
      }
    };
    class RangeFunction : public MacroFunction
    {
    public:
      RangeFunction() : MacroFunction("range") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        std::vector<std::string> secs=stringfx::split(stringfx::trim(in),',');
        if(!secs.size()) {cerr << "Bad input to function 'range'\n"; throw BadFunctionInputException();}
        double st=0,en;
        double step=1;
        if(secs.size()>1) {st=std::stod(secs[0]);en=std::stod(secs[1]);}
        else {st=0; en=std::stod(secs[0]);}
        if(secs.size()>2) step=std::stod(secs[2]);
        std::vector<double>* rexp=new std::vector<double>();
        if(st!=en)
        {
          if(st<en) {for(double i=st;i<en;i+=step) rexp->push_back(i);}
          else {for(double i=st;i>en;i+=step) rexp->push_back(i);}
        }
        return Variable("funcout_range","doublevector",rexp);
      }
    } RANGE;
    class InputFunction : public MacroFunction
    {
    public:
      InputFunction() : MacroFunction("input") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        std::string* sample=new std::string("");
        cout << in;
        cin >> *sample;
        return Variable("funcout_input","string",sample);
      }
    } INPUT;
    class SplitMoleculesFunction : public MacroFunction
    {
    public:
      SplitMoleculesFunction() : MacroFunction("splitmols") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        Molecule* mol=((Molecule*)std::stol(in));
        std::vector<Molecule*>* rexp=new std::vector<Molecule*>();
        int rno=-1;
        Molecule* curmol=nullptr;
        for(Atom* a : mol->getAtoms())
        {
          Atom* exp=new Atom(*a); exp->setResidueNumber(a->getResidueNumber()); exp->setResidue(a->getResidue());
          if(a->getResidueNumber()!=rno)
          {
            rno=a->getResidueNumber();
            curmol=new Molecule(exp);
            rexp->push_back(curmol);
            continue;
          }
          curmol->addAtom(exp);
        }
        return Variable("funcout_splitmols","molvector",rexp);
      }
    } SPLITMOLS;
    class MathEvaluateFunction : public MacroFunction
    {
    public:
      MathEvaluateFunction() : MacroFunction("matheval") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        //BODMAS
        int di=-1; char oc;
        std::string o1,o2;
        di=stringfx::lastIndexOf('[',in);
        while(di!=-1)
        {
          std::string repl=""; di++;
          while(in[di]!=']') repl+=in[di++];
          in=stringfx::replace(in,"["+repl+"]",evaluate(repl,m).getReplacementString());
          di=stringfx::lastIndexOf('[',in);
        }
        di=stringfx::indexOfAnyByPriority("^%/*+-",in);
        std::string er,reps;
        bool eonp=false;
        while(di!=std::string::npos)
        {
          o1=traceNumber(in,di-1,true); o2=traceNumber(in,di+1,false);
          reps=o1+std::string(1,in[di])+o2;
          switch(in[di])
          {
            case '^':
              er=to_string(::pow(std::stod(o1),std::stod(o2)));
              break;
            case '/':
              er=to_string(std::stod(o1)/std::stod(o2));
              break;
            case MODCHAR:
              er=to_string((int)std::stod(o1)%(int)std::stod(o2));
              break;
            case '*':
              er=to_string(std::stod(o1)*std::stod(o2));
              break;
            case '+':
              er=to_string(std::stod(o1)+std::stod(o2));
              break;
            case '-':
              if(!o1.length()) {eonp=true; break;}
              else er=to_string(std::stod(o1)-std::stod(o2));
              break;
          }
          if(eonp) break;
          in=stringfx::replace(in,reps,er);
          di=stringfx::indexOfAnyByPriority("^/*+-",in);
        }
        di=stringfx::indexOfAnyByPriority("<>=",in);
        while(di!=std::string::npos)
        {
          bool n1=true,n2=true; std::string rseg(1,in[di]);
          if(in[di+1]=='=')
          {
            rseg+='=';
            o2=traceNumber(in,di+2,false);
            if(stringfx::trim(o2)=="") {o2=traceString(in,di+2,false); n2=false;}
          }
          else
          {
            o2=traceNumber(in,di+1,false);
            if(stringfx::trim(o2)=="") {o2=traceString(in,di+1,false); n2=false;}
          }
          if(di>0 && in[di-1]=='!')
          {
            rseg='!'+rseg;
            o1=traceNumber(in,di-2,true);
            if(stringfx::trim(o1)=="") {o1=traceString(in,di-2,true); n1=false;}
          }
          else
          {
            o1=traceNumber(in,di-1,true);
            if(stringfx::trim(o1)=="") {o1=traceString(in,di-1,true); n1=false;}
          }
          reps=o1+rseg+o2;
          //cout<<"'" << reps << "'\n";
          switch(in[di])
          {
            case '>':
              if(in[di+1]=='=') er=to_string((int)(std::stod(o1)>=std::stod(o2)));
              else er=to_string((int)(std::stod(o1)>std::stod(o2)));
              break;
            case '<':
              if(in[di+1]=='=') er=to_string((int)(std::stod(o1)<=std::stod(o2)));
              else er=to_string((int)(std::stod(o1)<std::stod(o2)));
              break;
            case '=':
              if(di>0 && in[di-1]=='!')
              {
                if(n1 && n2) er=to_string((int)(std::stod(o1)!=std::stod(o2)));
                else if(!n1 && !n2) er=to_string((int)(stringfx::trim(o1)!=stringfx::trim(o2)));
                else er="1";
              }
              else
              {
                if(n1 && n2) er=to_string((int)(std::stod(o1)==std::stod(o2)));
                else if(!n1 && !n2) er=to_string((int)(stringfx::trim(o1)==stringfx::trim(o2)));
                else er="0";
              }
          }
          in=stringfx::replace(in,reps,er);
          di=stringfx::indexOfAnyByPriority("<>=",in);
        }
        di=stringfx::indexOfAnyByPriority("&|!",in);
        while(di!=std::string::npos)
        {
          o1=traceNumber(in,di-1,true); o2=traceNumber(in,di+1,false);
          reps=o1+std::string(1,in[di])+o2;
          switch(in[di])
          {
            case '!':
              er=((int)std::stod(o2))?"0":"1";
              break;
            case '&':
              er=((int)std::stod(o2) && (int)std::stod(o1))?"1":"0";
              break;
            case '|':
              er=((int)std::stod(o2) || (int)std::stod(o1))?"1":"0";
              break;
          }
          in=stringfx::replace(in,reps,er);
          di=stringfx::indexOfAnyByPriority("!&|",in);
        }
        return Variable("funcout_matheval","double",in);
      }
    } MATHEVAL;
    class DistanceSquare : public MacroFunction
    {
    public:
      DistanceSquare() : MacroFunction("distance2") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        std::string p1,p2;
        in=stringfx::trim(in);
        if(in[0]!='(')
        {
          cerr << "inbuiltfuncs::distance2 expected a pair of (double/int) vectors as input. Received: '"<<in<<"'\n";
          throw InvalidMathOperationException();
        }
        else
        {
          p1="";
          int ind=1;
          while(in[ind]!=')') p1+=in[ind++];
          int bloc=in.find("(",ind);
          if(bloc==std::string::npos)
          {
            cerr << "inbuiltfuncs::distance2 expected a pair of (double/int) vectors as input. Received: '"<<in<<"'\n";
            throw InvalidMathOperationException();
          }
          else
          {
            ind=bloc+1;
            p2="";
            while(in[ind]!=')') p2+=in[ind++];
          }
        }
        std::vector<string> ns1=stringfx::split(p1,','),ns2=stringfx::split(p2,',');
        if(ns1.size()!=ns2.size()) {cerr << "distance2: Vectors have unequal number of values\n"; throw InvalidMathOperationException();}
        double* r=new double(0);
        for(int i=0;i<ns1.size();i++) *r+=::pow(std::stod(ns1[i])-std::stod(ns2[i]),2);
        return Variable("funcout_distance2","double",r);
      }
    } DISTANCE2;
    class Distance : public MacroFunction
    {
    public:
      Distance() : MacroFunction("distance") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        double* num=(double*)(DISTANCE2.evaluate(in,m).getRawObjectPointer());
        *num=sqrt(*num);
        return Variable("funcout_distance","double",num);
      }
    } DISTANCE;
    class VectorBetweenAtoms : public ReferenceMacroFunction
    {
    public:
      VectorBetweenAtoms() : ReferenceMacroFunction("joiningvector") {}
      Variable evaluate(std::string in,const Macro& mc) override
      {
        std::vector<std::string> slcs=stringfx::split(in,',');
        Variable var=mc.getVar(stringfx::trim(slcs[0]));
        std::vector<double>* ret=new std::vector<double>();
        Molecule* m=((Molecule*)var.getRawObjectPointer());
        Eigen::Vector3d rv;
        try {rv=m->getAtoms()[std::stoi(slcs[1])]->bondVectorTo(m->getAtoms()[std::stoi(slcs[2])]);}
        catch(std::exception ex)
        {
          Atom *a1=nullptr,*a2=nullptr;
          slcs[1]=stringfx::trim(slcs[1]);
          slcs[2]=stringfx::trim(slcs[2]);
          for(Atom* ar : m->getAtoms())
          {
            if(!a1 && ar->toString()==slcs[1]) a1=ar;
            else if(!a2 && ar->toString()==slcs[2]) a2=ar;
            if(a1 && a2) break;
          }
          rv=a1->bondVectorTo(a2);
        }
        ret->push_back(rv(0)); ret->push_back(rv(1)); ret->push_back(rv(2));
        return Variable("funcout_joiningvector","doublevector",ret);
      }
    } JOININGVECTOR;
    class DotProductFunction : public ReferenceMacroFunction
    {
    public:
      DotProductFunction() : ReferenceMacroFunction("dot") {}
      Variable evaluate(std::string in,const Macro& mc) override
      {
        std::vector<double> r1,r2;
        auto vs=readTwoVectors(in,mc,"dot-product"); r1=get<0>(vs); r2=get<1>(vs);
        double* dp=new double(0);
        for(int i=0;i<r1.size();i++) *dp+=r1[i]*r2[i];
        return Variable("funcout_dot","double",dp);
      }
    } DOTPRODUCT;
    class CrossProductFunction : public ReferenceMacroFunction
    {
    public:
      CrossProductFunction() : ReferenceMacroFunction("cross") {}
      Variable evaluate(std::string in,const Macro& mc) override
      {
        std::vector<double> r1,r2;
        auto vs=readTwoVectors(in,mc,"cross-product",3); r1=get<0>(vs); r2=get<1>(vs);
        std::vector<double>* retv=new std::vector<double>(3);
        (*retv)[0]=r1[1]*r2[2]-r2[1]*r1[2];
        (*retv)[1]=r1[2]*r2[0]-r2[2]*r1[0];
        (*retv)[2]=r1[0]*r2[1]-r2[0]*r1[1];
        return Variable("funcout_cross","doublevector",retv);
      }
    } CROSSPRODUCT;
    class AngleFunction : public ReferenceMacroFunction
    {
    public:
      AngleFunction() : ReferenceMacroFunction("angleof") {}
      Variable evaluate(std::string in,const Macro& mc) override
      {
        std::vector<double> r1,r2;
        auto vs=readTwoVectors(in,mc,"angleof"); r1=get<0>(vs); r2=get<1>(vs);
        double *dp=new double(0),m1=0,m2=0;
        for(int i=0;i<r1.size();i++) {*dp+=r1[i]*r2[i]; m1+=r1[i]*r1[i]; m2+=r2[i]*r2[i];}
        *dp=acos(*dp/(sqrt(m1)*sqrt(m2)));
        return Variable("funcout_angle","double",dp);
      }
    } ANGLEBETWEEN;
    class DihedralFunction : public ReferenceMacroFunction
    {
    public:
      DihedralFunction() : ReferenceMacroFunction("dihedralof") {}
      Variable evaluate(std::string in,const Macro& mc) override
      {
        std::vector<double> r1,r2,r3;
        auto vs=readManyVectors(in,mc,"dihedralof"); r1=vs[0]; r2=vs[1]; r3=vs[2];
        double *dp=new double(0),m1=0,m2=0;
        std::vector<double> n1(3),n2(3);
        n1[0]=r1[1]*r2[2]-r2[1]*r1[2];
        n1[1]=r1[2]*r2[0]-r2[2]*r1[0];
        n1[2]=r1[0]*r2[1]-r2[0]*r1[1];
        n2[0]=r2[1]*r3[2]-r3[1]*r2[2];
        n2[1]=r2[2]*r3[0]-r3[2]*r2[0];
        n2[2]=r2[0]*r3[1]-r3[0]*r2[1];
        for(int i=0;i<r1.size();i++) {*dp+=n1[i]*n2[i]; m1+=n1[i]*n1[i]; m2+=n2[i]*n2[i];}
        *dp=acos(*dp/(sqrt(m1)*sqrt(m2)));
        return Variable("funcout_dihedral","double",dp);
      }
    } DIHEDRALBETWEEN;
    class ArrayAppend : public ReferenceMacroFunction
    {
    public:
      ArrayAppend() : ReferenceMacroFunction("append") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        std::vector<string> ins=stringfx::split(in,',',false,'(',')');
        Variable v=m.getVar(ins[0]);
        if(m.hasVar(ins[1])) ((std::vector<Variable>*)v.getRawObjectPointer())->push_back(m.getVar(ins[1]));
        else
        {
          int spc=ins[1].find(' '),brk=ins[1].find('(');
          if(spc!=std::string::npos && (brk==std::string::npos || brk>spc))
          {
            std::stringstream ss(ins[1]);
            std::string tp; ss >> tp;
            std::string rest=stringfx::drain(ss);
            ((std::vector<Variable>*)v.getRawObjectPointer())->push_back(Variable("array_ref",tp,rest,m.getMacroFF()));
          }
          else ((std::vector<Variable>*)v.getRawObjectPointer())->push_back(Variable("array_ref",v.getTypeString(),ins[1],m.getMacroFF()));
        }
        return Variable("rv","string","pass");
      }
    } APPENDARRAY;
    class MultiMolLoad : public MacroFunction
    {
    public:
      MultiMolLoad() : MacroFunction("loadtraj") {} //Implicitly GRO format
      Variable evaluate(std::string in,const Macro& m) override
      {
        std::ifstream ldf; ldf.open(stringfx::trim(in));
        std::vector<Molecule*> *dt=new std::vector<Molecule*>(),ld=molfx::loadManyGROMolecules(ldf,*(m.getMacroFF()));
        for(int i=0;i<ld.size();i++) dt->push_back(ld[i]);
        ldf.close();
        return Variable("funcout_loadtraj","molvector",dt);
      }
    } LOADTRAJ;
    /*class CondEvaluationFunction : public MacroFunction
    {
    public:
      CondEvaluationFunction() : MacroFunction("condeval") {}
      Variable evaluate(std::string in,const Macro& m) override
      {
        int di=-1; char oc;
        std::string o1,o2;
        di=stringfx::lastIndexOf('[',in);
        while(di!=-1)
        {
          std::string repl=""; di++;
          while(in[di]!=']') repl+=in[di++];
          in=stringfx::replace(in,"["+repl+"]",evaluate(repl,m).getReplacementString());
          di=stringfx::lastIndexOf('[',in);
        }
      }
    } CONDEVAL;*/
    static void init(Macro& m) {m.addFunction(&RANGE); m.addFunction(&INPUT); m.addFunction(&LOADTRAJ); m.addFunction(&SPLITMOLS); m.addFunction(&MATHEVAL); m.addFunction(&DISTANCE2); m.addFunction(&DISTANCE); m.addFunction(&APPENDARRAY); m.addFunction(&JOININGVECTOR); m.addFunction(&DOTPRODUCT); m.addFunction(&CROSSPRODUCT); m.addFunction(&ANGLEBETWEEN); m.addFunction(&DIHEDRALBETWEEN);}
  }
  class Prompt
  {
    Macro mymacro;
    std::string prompt;
  public:
    Prompt(const std::string& pr="") {prompt=pr; inbuiltfuncs::init(mymacro); launchPrompt();}

    void launchPrompt()
    {
      std::string iln;
      bool succ=true;
      while(succ)
      {
        if(cin.eof()) break;
        cout << prompt<<">\t";
        getline(cin,iln);
        succ=mymacro.addLineandRun(iln);
      }
    }
  };
}
//Completing the variable class
void macros::Variable::assign(const std::string& ls)
{
  std::vector<int>* iv=nullptr;
  std::vector<double>* dv=nullptr;
  std::vector<std::string>* sv=nullptr;
  std::vector<Molecule*>* mv=nullptr;
  std::vector<macros::Variable>* vv=nullptr;
  std::string lop=""; char ch; bool begin=false;
  int stg=-1,K=0;
  switch(type)
  {
    case 0:
      if(data) delete ((int*)data);
      if(!ls.length()) data=nullptr;
      else data=new int(std::stod(ls));
      break;
    case 1:
      if(data) delete ((double*)data);
      if(!ls.length()) data=nullptr;
      else data=new double(std::stod(ls));
      break;
    case 2:
      if(data) delete ((std::string*)data);
      if(!ls.length()) data=nullptr;
      else data=new std::string(ls);
      break;
    case 3:
      if(data) delete ((Atom*)data);
      if(!ls.length()) data=nullptr;
      else
      {
        try {data=(Atom*)std::stol(ls);}
        catch(std::invalid_argument e) {data=new Atom(stringfx::split(ls,' ')[0],0,0,0);}
      }
      break;
    case 4:
      if(data) delete ((Molecule*)data);
      if(!ls.length()) data=nullptr;
      else
      {
        try {data=(Molecule*)std::stol(ls);}
        catch(std::exception e) {data=new Molecule(stringfx::split(ls,' ')[0],*myff,stringfx::split(ls,' ')[1]);}
      }
      break;
    case 5:
      if(data) delete (std::vector<int>*)data;
      if(!ls.length()) {data=nullptr; break;}
      iv=new std::vector<int>();
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {begin=true; continue;}
        if(ch==')') {begin=false; break;}
        if(begin && (ch>='0' && ch<='9') || ch=='-') lop+=ch;
        if(ch==',') {iv->push_back((int)std::stod(lop)); lop="";}
      }
      if(lop.length()) iv->push_back((int)std::stod(lop));
      data=iv;
      break;
    case 6:
      if(data) delete (std::vector<double>*)data;
      if(!ls.length()) {data=nullptr; break;}
      dv=new std::vector<double>();
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {begin=true; continue;}
        if(ch==')') {begin=false; break;}
        if(begin && ((ch>='0' && ch<='9') || ch=='-' || ch=='.')) lop+=ch;
        if(ch==',') {dv->push_back(std::stod(lop)); lop="";}
      }
      if(lop.length()) dv->push_back(std::stod(lop));
      data=dv;
      break;
    case 7:
      if(data) delete (std::vector<std::string>*)data;
      if(!ls.length()) {data=nullptr; break;}
      sv=new std::vector<std::string>();
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {begin=true; continue;}
        if(ch==')') {begin=false; break;}
        if(begin) lop+=ch;
        if(ch==',') {sv->push_back(lop); lop="";}
      }
      if(lop.length()) sv->push_back(lop);
      data=sv;
      break;
    case 8:
      if(data) delete (std::vector<Molecule*>*)data;
      if(!ls.length()) {data=nullptr; break;}
      mv=new std::vector<Molecule*>();
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {stg++; if(!stg) {begin=true; continue;}}
        if(ch==')') {stg--; if(stg==-1) {begin=false; break;}}
        if(ch==',' && !stg) {mv->push_back((Molecule*)std::stol(lop)); lop=""; continue;}
        if(begin) lop+=ch;
      }
      if(lop.length()) mv->push_back((Molecule*)std::stol(lop));
      data=mv;
      break;
    case 9:
      if(data) delete (std::vector<macros::Variable>*)data;
      if(!ls.length()) {data=nullptr; break;}
      vv=new std::vector<macros::Variable>();
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {stg++; if(!stg) {begin=true; continue;}}
        if(ch==')') {stg--; if(stg==-1) {begin=false; break;}}
        if(ch==',' && !stg) {vv->push_back(macros::Variable("arrayindex",typestr,lop,myff)); lop=""; continue;}
        if(begin) lop+=ch;
      }
      if(lop.length()) vv->push_back(macros::Variable("arrayindex",typestr,lop,myff));
      data=vv;
      break;
    case 10:
      if(data) delete (std::vector<macros::Variable>*)data;
      if(!ls.length()) {data=nullptr; break;}
      vv=new std::vector<macros::Variable>();
      K=0;
      for(int i=0;i<ls.length();i++)
      {
        ch=ls[i];
        if(ch=='(') {stg++; if(!stg) {begin=true; continue;}}
        if(ch==')') {stg--; if(stg==-1) {begin=false; break;}}
        if(ch==',' && !stg) {vv->push_back(macros::Variable("arrayindex",myStructTypes[K++],lop,myff)); lop=""; continue;}
        if(begin) lop+=ch;
      }
      if(lop.length()) vv->push_back(macros::Variable("arrayindex",myStructTypes[K],lop,myff));
      data=vv;
      break;
    default:
      cerr << "Could not assign '"<<ls<<"' to an unknown typeID: "<<type<<"\n";
      throw BadTypeException();
  }
}
inline std::string macros::Variable::getReplacementString() const
{
  if(!data) return "";
  std::string ret="";
  int vl;
  switch(type)
  {
    case 0:
      return to_string(*(int*)data);
    case 1:
      return to_string(*(double*)data);
    case 2:
      return (*(std::string*)data);
    case 3:
    case 4:
      return to_string((long)data);
    case 5:
      ret="(";
      vl=((std::vector<int>*)data)->size();
      if(vl) ret+=to_string(((std::vector<int>*)data)->operator[](0));
      for(int i=1;i<vl;i++) ret+=(","+to_string(((std::vector<int>*)data)->operator[](i)));
      ret+=')';
      return ret;
    case 6:
      ret="(";
      vl=((std::vector<double>*)data)->size();
      if(vl) ret+=to_string(((std::vector<double>*)data)->operator[](0));
      for(int i=1;i<vl;i++) ret+=(","+to_string(((std::vector<double>*)data)->operator[](i)));
      ret+=')';
      return ret;
    case 7:
      ret="(";
      vl=((std::vector<std::string>*)data)->size();
      if(vl) ret+=((std::vector<std::string>*)data)->operator[](0);
      for(int i=1;i<vl;i++) ret+=(","+((std::vector<std::string>*)data)->operator[](i));
      ret+=')';
      return ret;
    case 8:
      ret="(";
      vl=((std::vector<Molecule*>*)data)->size();
      if(vl) ret+=to_string((long)((std::vector<Molecule*>*)data)->operator[](0));
      for(int i=1;i<vl;i++) ret+=(","+to_string((long)((std::vector<Molecule*>*)data)->operator[](i)));
      ret+=')';
      return ret;
    case 9:
    case 10:
      ret="(";
      vl=((std::vector<macros::Variable>*)data)->size();
      if(vl) ret+=((std::vector<macros::Variable>*)data)->operator[](0).getReplacementString();
      for(int i=1;i<vl;i++) ret+=(","+((std::vector<macros::Variable>*)data)->operator[](i).getReplacementString());
      ret+=')';
      return ret;
    default:
      throw macros::BadTypeException();
  }
}
void macros::Variable::clearData()
{
  switch(type)
  {
    case 0:
      delete (int*) data;
      break;
    case 1:
      delete (double*) data;
      break;
    case 2:
      delete (std::string*) data;
      break;
    case 3:
      delete (Atom*) data;
      break;
    case 4:
      delete (Molecule*) data;
      break;
    case 5:
      delete (std::vector<int>*) data;
      break;
    case 6:
      delete (std::vector<double>*) data;
      break;
    case 7:
      delete (std::vector<std::string>*) data;
      break;
    case 8:
      delete (std::vector<Molecule*>*) data;
      break;
    case 9:
    case 10:
      delete (std::vector<macros::Variable>*) data;
      break;
  }
}

//Completing Macro
void macros::Macro::callMacro(const std::string& cm)
{
  Macro mn(cm,*this);
  for(MacroFunction* mf : funcs) mn.addFunction(mf);
  mn.execMacro();
}
#endif
