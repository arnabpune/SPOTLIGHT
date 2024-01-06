#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "../graph/Atom.hpp"
//#include "../support/commons.h"

using namespace std;

//Exceptions
/**This exception is throws when any of the required parameters are not available*/
class DataNotAvailableException: public exception {};
/**This exception is thrown when a function called on a bond object does not contain an atom it is expecting. See BondData*/
class AtomNotInBondException :  public exception {};
/**This exception is thrown when no valid atom-type is available to bond to an atom in the bond data*/
class NoBondsAvailableException :  public exception {};
/**This exception is thrown when the rule criterion is extremely strict. This happens when the number of valencies for an atom-type are less than the minimum number required by the rules*/
class RulesExceedValenciesException : public exception {}; //If the number of rules for an atom type are more than its valences
/**This exception is thrown when no atom-type can bond to a given atom due to rule limitations.<br/> This is different from NoBondsAvailableException*/
class RuleConstraintsPreventBondsException : public exception {}; //If the rules are so restrictive that no atom bond is allowed
/**This exception is thrown when a requested category is not present in the force-field*/
class NoSuchCategoryException : public exception {}; //Self-explanatory
/**This exception is thrown when a residue is not in the list while loading a Protein*/
class ResidueNotFoundException : public exception {}; //Self-explanatory
/**This exception is thrown when a residue being loaded contains an atom type that is not registered in that residue*/
class AtomNotInResidueException :  public exception {};
//Some typedefs
typedef std::vector<double> vector_d;
typedef std::pair<double,double> pair_dd;
typedef std::pair<std::string,std::string> pair_ss;
typedef std::vector<std::string> vector_s;
//Common simple functions
template<class T> inline bool matches(const std::pair<T,T>& p1,const std::pair<T,T>& p2) {return (((get<0>(p1))==(get<0>(p2)) && (get<1>(p1))==(get<1>(p2))) || ((get<0>(p1))==(get<1>(p2)) && (get<1>(p1))==(get<0>(p2))));}
template<class T> bool matches(const std::vector<T>& l1,const std::vector<T>& l2)
{
	if(l1.size()!=l2.size()) return false;
	bool m1=true,m2=true;
	for(int i=0;i<l1.size();i++)
	{
		if(l1[i]!=l2[i]) m1=false;
		if(l1[i]!=l2[l1.size()-1-i]) m2=false;
		if(!m1 && !m2) return false;
	}
	return (m1||m2);
}
//Useful constants
#ifdef STATICDATA
	#ifdef AMINOACIDDATA
		static const std::string ff_folder="DNV_ROOT/data/aminoacids/itps"; //Change accordingly
	#else
		static const std::string ff_folder="DNV_ROOT/data/itps"; //Change accordingly
	#endif
#else
	static const std::string ff_folder="./itps"; //Change accordingly
#endif
static const std::string LengthDataFile=ff_folder+"/bondtypes.itp",AngleDataFile=ff_folder+"/angletypes.itp",DihedralDataFile=ff_folder+"/dihedraltypes.itp",ImproperDataFile=ff_folder+"/impropertypes.itp";
static const std::string OneFourFile=ff_folder+"/pairtypes.itp",OtherInterFile=ff_folder+"/nonbonded.calc.itp";
static const double THETACUT=(7*PI)/180.0; /**<Default angle fluctuation tolerance (7°)*/
static const double DIHEDTHETACUT=(15*PI)/180.0; /**<Default dihedral fluctuation tolerance (15°)*/

//Energy functions
inline static double dihedralEnergyFunction(const std::vector<double>& par,double ph) {return par[1]*(1+cos(par[2]*ph-par[0])); }
inline static double angleEnergyFunction(const std::vector<double>& par,double th,double ub) {return (((par[0]-th)*(par[0]-th)*par[1]) +((par[2]-ub)*(par[2]-ub)*par[3]))/2;}

class ForceField;
/**\class Residue forcefield.hpp "forcefield/forcefield.hpp"
	This is the residue class. It loads in residues which are building blocks for loading standard molecules (especially proteins).<br/> This residue data is available in the RTP files (eg. protein.rtp, dna.rtp) and these files are available with DeNovo even though they were extracted from GROMACS
	Residue names are recognized while loading from PDB or GRO formats and the atom-names are mapped accordingly.<br/>It stores the different atom-names (called aliases) and the corresponding atom-type names as well.
*/
class Residue
{
	std::vector<std::string> ids,actual;
	std::vector<double> charges;
	std::string name;
public:
	/**@brief Construct a residue by providing DeNovo with the atom-type definition lines from the RTP file*/
	Residue(const std::string& rnm,const std::vector<std::string>& lines)
	{
		name=rnm;
		ids=std::vector<std::string>(); actual=std::vector<std::string>();
		std::string nm,rl,w1,w2;
		double q;
		for(const std::string& s : lines)
		{
			stringstream ss(s);
			ss >> nm >> rl >> q;
			ids.push_back(nm); actual.push_back(rl); charges.push_back(q);
		}
	}

	/**@brief Get the number of atoms defined in this Residue*/
	inline int getAtomCount() const {return ids.size();}
	/**@brief Get the list of aliases (unique atom labels given for the residue based on its position)*/
	inline const std::vector<std::string>& getAliases() const {return ids;}
	/**@brief Get the list of actual atom-type names in the same order as they were loaded into the residue from the RTP file*/
	inline const std::vector<std::string>& getActuals() const {return actual;}
	/**@brief Get the list of actual charges (reassigned) in the same order as they were loaded into the residue from the RTP file*/
	inline const std::vector<double>& getAssignedCharges() const {return charges;}
	/**@brief Get the 'x'th alias*/
	inline std::string getAlias(int x) const {return ids[x];}
	/**@brief Get the 'x'th actual atom-type name*/
	inline std::string getActual(int x) const {return actual[x];}
	/**@brief Get the name of the Residue*/
	inline const std::string& getName() const {return name;}

	/**@brief Get the actual (force-field) atom-type name for a given alias*/
	std::string getActual(const std::string& alias) const
	{
		for(int i=0;i<ids.size();i++) {if(ids[i]==alias) return actual[i];}
		cout << alias << " not in "<<name<<"\n";
		throw AtomNotInResidueException();
	}
	/**@brief Get the actual (force-field) atom-type along with reassigned charge (if applicable)*/
	std::pair<std::string,double> getActualWithCharge(const std::string& alias) const
	{
		for(int i=0;i<ids.size();i++) {if(ids[i]==alias) return make_pair(actual[i],charges[i]);}
		cerr << alias << " not in "<<name<<"\n";
		throw AtomNotInResidueException();
	}
};
//Bond
/**\class BondData forcefield.hpp "forcefield/forcefield.hpp"
	This is the bond data class. It stores the bond parameter data (atom-type pairs that form the bond, the bond length and force constant)<br/>
	See the documentation (PDF) for more details on bond data format
*/
class BondData
{
public:
	std::string atom1; /**<The first atom-type which forms the bond (according to file - actual order doesn't matter)*/
	std::string atom2; /**<The second atom-type which forms the bond (according to file - actual order doesn't matter)*/
	double b0; /**<Equlibrium bond length*/
	double kb; /**<The force-constant*/

	/**@brief The standard constructor for BondData which takes in two atom-type names and the two bond-parameters*/
	BondData(std::string at1,std::string at2,double p1,double p2) {atom1=at1; atom2=at2; b0=p1; kb=p2;}

	/**@brief Get the pair of atoms (as std::pair<std::string,std::string>) forming this bond*/
	std::pair<std::string,std::string> getAtoms() const {return make_pair(atom1,atom2);}
	/**@brief Get the bond parameters (as a pair of doubles - std::pair<double,double>)*/
	std::pair<double,double> getParameters() const {return make_pair(b0,kb);}
	/**@brief Checks if the given atom-type corresponds to either of the two atom-types which contribute to this bond*/
	inline bool contains(const std::string& s) const {return (atom1==s || atom2==s);}
	/**@brief Checks if the BondData object is composed of the same two atom-types as the two provided*/
	inline bool satisfies(const std::string& s1,const std::string& s2) const {return (atom1==s1 && atom2==s2) || (atom1==s2 && atom2==s1);}
	/**@brief Get the second atom in the pair forming this bond (given one of the two)*/
	const std::string& getSecond(const std::string& s) const
	{
		if(atom1==s) return atom2;
		if(atom2==s) return atom1;
		throw AtomNotInBondException();
	}
	/**@brief Get the second atom in the pair forming this bond (given one of the two)*/
	inline Atom* getSecond(Atom* a,const ForceField& ff) const;
	/**@brief Get the equilibrium length parameter for this bond*/
	inline double length() const {return b0;}
	/**@brief Get the equilibrium length parameter for this bond*/
	inline double getLength() const {return b0;}
	/**@brief Get a printable string representation of this BondData object*/
	inline std::string toString() const {return atom1+" "+atom2+"\t("+to_string(b0)+","+to_string(kb)+")";}
	/**@brief Replace one atom in this BondData object by another
		 @param[in] orig: The original atom which belongs to this object which must be replaced
		 @param[in] rep: The replacement atom-type name
	*/
	void replace(const std::string& orig,const std::string& rep)
	{
		if(atom1==orig) atom1=rep;
		else if(atom2==orig) atom2=rep;
		else throw AtomNotInBondException();
	}
};
//End Bond
namespace chemtools
{
	/**\class Rule forcefield.hpp "forcefield/forcefield.hpp"
		 The rule class stores condition data which is used during generation to ensure that the various atom-types satisfy their definitions.
		 The exact functioning of the rules is given in the documentation (PDF)<br/>
		 Rules are defined on specific atom-types (referred to as the central atom) and satisfying these rules is supposed to ensure that the atom-type will mimic all properties as in the forcefield where it was parameterized.<br/>
		 Some extra parameters (which are now not recommended to be used) are not mentioned in documentation: They are angle and dihedral forcing. The rule would force a particular angle or dihedral to be used (rather than the standard forcefield parameters). This feature is deprecated since DeNovo θ.
	*/
	class Rule
  {
    std::string centre;
    std::vector<std::string> group;
    int min=0,max=10;
		double angv,dihv,dihm;
	public:
		bool angp=false; /**<\deprecated Has a preferred angle?*/
		bool dihp=false; /**<\deprecated Has a preferred dihedral?*/

  public:
		/**@brief The blank constructor. Only to be used as a placeholder*/
    Rule() : Rule("DUM","DUM") {}
		/**@brief The Rule constructor. It creates the Rule object taking input in a standard format
			 @details It can take only one rule at a time and cannot directly interpret Category names (See buildRules() which can take the raw line from the file). A single rule ends with a semicolon ';' or a newline '\n' whichever comes first<br/>
			 A sample rule (singular) looks like this
			 @code{.txt}
			 	H* 1 1
				OH1|FA|CF3 1 1
			 @endcode
			 @param[in] src: The source atom for which the rule is meant to apply
			 @param[in] rule: The single rule (in string format)
			 @param[in] ff: The ForceField from which various atom-type names used in the rules are referred.
		*/
    Rule(std::string src,std::string rule,const ForceField* ff=nullptr);

		/**@brief Checks if rules permit adding in a new atom
			 @details This functions checks to see if rules on the atom-types described by 's' and 'n' allow 'n' (a new atom with no bonds upto now) to bond to a source atom 's'.<br/>
			 This function returns negative (false) if bonding these two atoms would violate a rule constraint (eg. the maximum limit of bonds allowed for 's' or 'n' with the other).
			 @param[in] s: The source atom which is already present in the molecule to which the new atom is to be attached
			 @param[in] n: The new atom which is being trialled. It has no bonds upto now.
			 @param[in] m: The molecule from which 's' comes. This will be used to find all the atoms bonded to 's'
			 <br/>See also: satisfies(Atom*,Molecule*,const ForceField*)
		*/
    inline bool satisfies(Atom* s,Atom* n,Molecule* m) const;
		/**@brief Alternate method for satisfies(Atom*,Atom*,Molecule*) which uses the bonded atom list rather than the Molecule to calculate it*/
    bool satisfies(Atom* s,Atom* n,const std::vector<Atom*>& b) const;
		/**@brief Tests if an atom satisfies the requirements of this Rule in it's current state (based on what it is bonded to)*/
		inline bool satisfies(Atom* s,Molecule* m,const ForceField* ff=nullptr) const;
		/**@brief Alternate method for satisfies(Atom*,Molecule*,const ForceField*) which uses the bonded atom list rather than the Molecule to calculate it*/
		bool satisfies(Atom* s,const std::vector<Atom*>& bnd,const ForceField* ff=nullptr) const;
		/**@brief \deprecated Same as satisfies(Atom*,Atom*,Molecule*) but does the calculation assuming that the extra bond between the source and new atom is also present*/
    inline bool satisfiesStrong(Atom* s,Atom* n,Molecule * m) const;
		/**@brief \deprecated An alternate function for satisfiesStrong(Atom*,Atom*,Molecule*)*/
    bool satisfiesStrong(Atom* s,Atom* n,const std::vector<Atom*>& b) const;

		/**@brief Get the atom-type name on which this Rule is defined (an alternative for a "const" Rule object)*/
		inline const std::string& getCentre() const {return centre;}
		/**@brief Get the atom-type name on which this Rule is defined*/
    std::string getCentre() {return centre;}
		/**@brief Get the group of atoms which this rule mentions.  (an alternative for a "const" Rule object)*/
    inline const std::vector<std::string>& getRuleGroup() const {return group;}
		/**@brief Get the group of atoms which this rule mentions.
			 @details This does not include the centre atom (see getCentre())<br/>
			 This will return all the atom-types to which the central atom can bond which are mentioned in this rule.
		*/
    std::vector<std::string>& getRuleGroup() {return group;}
		/**@brief Get the maximum number of bonds allowed by this Rule with the atom-types listed in it.*/
    inline int getMax() const {return max;}
		/**@brief Get the minimum number of bonds required by this Rule with the atom-types listed in it.*/
    inline int getMin() const {return min;}
		inline double getAngle() const {return angv;}
		inline double getDihedral() const {return dihv;}
		/**@brief Get the angle specified by this rule
			 @details Older versions of DeNovo allowed specifying an angle to override the force-field parameters when bonding to this atom as the center.<br/>
			 \deprecated It is deprecated since DeNovo θ (at-least). Support is not guaranteed for previous versions and absent in newer versions.
		*/
		inline std::vector<double> predictAngle(const ForceField& ff,Atom* a1,Atom* a2,Atom* a3);
		/**@brief Get the dihedral specified by this rule
			 @details Older versions of DeNovo allowed specifying an dihedral value to override the force-field parameters when bonding to this atom as the center.<br/>
			 \deprecated It is deprecated since DeNovo θ (at-least). Support is not guaranteed for previous versions and absent in newer versions.
		*/
		inline std::vector<double> predictDihedral(const ForceField& ff,Atom* a1,Atom* a2,Atom* a3,Atom* a4);
		/**@brief Choose a random atom-type from the list stored in this Rule.
			 @details This function returns a random atom that will (according to this Rule) satisfy the requirements of the central atom.<br/>When multiple atoms are specified using the '|' separator or a Category, it is assumed that any of these atoms can be used to satisfy the rule and selection of each is equally likely<br/>Atom-type reweighting is likely to work at this stage as well, but hasn't been tested so far
		*/
		inline Atom* suggest(Atom* c,const ForceField* ff) const;
  };

	/**@brief Generate a set of Rule objects by parsing a line from the definitions file (See PDF for definition format)
		 @details This function takes the raw line as input and splits the line into multiple rules (if present), and returns the entire list of rules.
		 <br/>Sample input lines look like this:
		 @code{.txt}
		 	H* 1 1;BRAL 0 1;CLAL 0 1
			[aromatic] 2;OH1|FA|CF3 1 1
		 @endcode
	*/
  static std::vector<Rule> buildRules(const std::string& def,const ForceField* ff=nullptr);
	/**@brief Load all rules from a file
		 @details Load all the rules atomtype-wise from the input file. buildRules() is used to build the rules.<br/>
		 See the documentation PDF for the atom-type definition file format
		 @returns a list with pairs of atom-types and the list of rules it is expected to satisfy (std::vector<std::pair<std::string,std::vector<Rule>>>)
	*/
  static std::vector<std::pair<std::string,std::vector<Rule>>> getRuleset(const std::string& file,const ForceField* ff=nullptr);
}
/** \class Category forcefield.hpp "forcefield/forcefield.hpp"
	 This class is used to store lists of atom-type names grouped under a category name. A single atom-type can be grouped under multiple categories.<br/>The format and working are briefly explained the in algorithm PDF.
	 <br/>An example category line looks like this:
	 @code{.txt}
	 	halogen FA|F1|F2|F3|CLAL|CLAR|BRAR|BRAL|I
	 @endcode
	 This is the "halogen" Category.
*/
class Category
{
	std::string targa;
	std::string reps;
	std::vector<std::string> myGroup;
public:
	/**@brief The standard (and simple) Category class constructor
		 @details It takes in a category name, and the '|'-separated list of atom names to construct this Category object
		 @param[in] s1: The name of the Category
		 @param[in] rl: The list of atom-type names contained in this Category (with '|' as the separator)
		 @param[in] ff: The forcefield from which the atom-type names of this Category come.
	*/
	Category(const std::string& s1,const std::string& rl,const ForceField* ff=nullptr);

	/**@brief Get the name of this Category*/
	inline const std::string& getCategoryName() const {return targa;}
	/**@brief Describe this Category in a string. ie. get a representative string (useful to print out)*/
	inline const std::string& toString() const {return reps;}
	/**Get a list of strings containing the atom-type names for the atom-types grouped under this category*/
	inline const std::vector<std::string>& getGroup() const {return myGroup;}
	/**@brief Checks whether a particular atom type is grouped under this Category*/
	inline bool contains(const Atom* a) const {return contains(a->toString());}
	/**@brief Checks whether a particular atom type is grouped under this Category (uses atom-type name instead of Atom object pointer)*/
	inline bool contains(const std::string& str) const {return ::contains(myGroup,str);}
	/**@brief Get the number of elements in this category (including repeats if any)*/
	inline int getSize() const {return myGroup.size();}
};
/**\class ForceField forcefield.hpp "forcefield/forcefield.hpp"
	 This is the ForceField class which contains all the force-field releated parameters, such as bond-lengths, angles, dihedrals, etc. It also contains the atom-type definition rules, and user-defined categories.<br/>
	 This class provides all the functions to quickly access this data by storing it in an organized format.<br/>
	 Various forcefields can be loaded into DeNovo if required parameters are available. Special cases where this is not true include the "empty" forcefield and the IncompleteForceField<br/>
	 <br/>See the documentation PDF for details about various forcefield parameters and the way they are stored and loaded, and the different forcefields like the "empty" forcefield and the "CHARMM-27" forcefield
*/
class ForceField
{
protected:
  std::vector<Atom*> atom_types; /**<The list that stores one Atom object copy for each atom-type. These objects are cloned (see Atom::Atom(const Atom&)) when a new Atom of a particular atom-type name is retrieved from this ForceField*/
	std::vector<Category> categories; /**<Stores the list of categories loaded into this ForceField (see: Category)*/
	std::vector<std::pair<std::string,std::vector<chemtools::Rule>>> rules; /**<Maps each atom type to a set of rules (which are the atom-type definition rules). See chemtools::Rule and atom-type definition rules in the documentation (PDF)*/
	std::vector<Residue> residues; /**<Residue objects can also be stored when loading residues from an RTP file (see loadResidues())*/
	std::vector<std::pair<char,double>> electronegativities; /**<Still under development - intended to be used for Gasteiger Charge calculation*/
	std::vector<std::pair<char,int>> defaultconnectivities; /**<Default connectivity data to be used (if no valencies are loaded - such as for an incomplete forcefield). Especially used for calculating Bond orders*/
public:
  std::vector<BondData> length_parameters; /**<Bond data stored in BondData objects. They contain a pair of atom-types which contribute to forming the bond, and the bond parameters (equilibrium length and force-constant) */
  std::vector<std::pair<vector_s,vector_d>> angle_parameters; /**<Angle parameter data. It maps string vectors (list of atoms contributing to the angle formation) to angle parameters (including Urey-Bradley parameters - See documentation PDF for details)*/
  std::vector<std::pair<vector_s,vector_d>> dihedral_parameters; /**<Dihedral parameter data. It maps string vectors (list of atom names for the atoms contributing to the dihedral formation in proper order) to dihedral parameters (equilibrium value,force-constant,multiplicity)*/
  std::vector<std::pair<vector_s,pair_dd>> improper_parameters; /**<Similar to dihedral_parameters. It maps atom lists to improper dihedral parameters*/
  std::vector<std::pair<std::string,std::vector<BondData>>> atombondmap; /**<This list maps each atom-type to a list of bonds (stored as BondData objects) which are bonds it is allowed to make.*/

  std::vector<std::pair<pair_ss,pair_dd>> one_four_interactions; /**<Stores 1-4 interaction parameters. Not used in later versions of DeNovo*/
  std::vector<std::pair<std::string,vector_d>> other_interactions; /**<Stores L-J interaction parameters. Not used in later versions of DeNovo. These parameters are taken directly from the atom-types*/
  //std::vector<std::pair<std::string,double>> one_four_charges;

public:
	/**@brief Blank constructor for ForceField (does not load any force-field parameters)
		 @details An empty constructor which loads no parameters except default values for standard elements (such as 4 for carbon, 3 for nitrogen, etc.).
		 It is to be called by subclass which loads the actual atom-types (such as IncompleteForceField). Ensure that any class that is derived from ForceField and uses this constructor defines the atom-types (with proper force-field energy parameters - σ,ε,q)
	*/
	ForceField()
	{
		if(!defaultconnectivities.size())
		{
			defaultconnectivities.push_back(make_pair('C',4));
			defaultconnectivities.push_back(make_pair('N',3));
			defaultconnectivities.push_back(make_pair('O',2));
			defaultconnectivities.push_back(make_pair('H',1));
			defaultconnectivities.push_back(make_pair('F',1));
			defaultconnectivities.push_back(make_pair('I',1));
			defaultconnectivities.push_back(make_pair('P',5));
			defaultconnectivities.push_back(make_pair('S',6));
		}
	}
public:
	/**@brief Load the forcefield from the standard ffin format (see PDF documentation for format details)
		 @details This function completely loads <b>all the required</b> parameters for the atom-types by reading data from the file. If this wasn't your requirement, or if data is not available, check the IncompleteForceField class.<br/>
		 This function uses the ForceField::loadAdditionalAtomTypes() functionality.
	*/
  ForceField(const std::string& ff_file) :ForceField()
  {
    atom_types=std::vector<Atom*>();
		cout << "Loading Forcefield from: "<<ff_file<<"\n";
    loadAdditionalAtomTypes(ff_file);
    cout <<"ForceField loaded! "<<atom_types.size()<<" atom types available and have all necessary parameters.\n";
  }
	/**@brief Load the forcefield from the standard ffin format, and categories from a categories file.  (see PDF documentation for format details)*/
	ForceField(const std::string& ff_file,const std::string& catf) : ForceField(ff_file)  {loadCategories(catf);}
	ForceField(const std::string& ff_file,const std::string& catf,const std::string& chgfile) : ForceField(ff_file,catf)  {loadChargeGroups(chgfile);}

	inline void loadResidues(const std::string& str); //Documented below
	/**@brief Get a list with all the residues loaded into this Forcefield*/
	inline const std::vector<Residue>& getResidues() const {return residues;}
	/**@brief Get a Residue from those loaded into this ForceField by name*/
	const Residue& getResidue(const std::string& rname) const
	{
		for(const Residue& r : residues) {if(r.getName()==rname) return r;}
		throw ResidueNotFoundException();
	}
	/**@brief Get the list of all Categories loaded into this ForceField*/
	inline const auto& /*std::vector<Category>&*/ getCategories() const {return categories;}
	/**@brief Get all the atom-types (as Atom object pointers) loaded into this ForceField*/
	inline const std::vector<Atom*>& getAtomTypes() const {return atom_types;}
	/**@brief Get the defaualt connectivity (not from forcefield data) for a given element*/
	inline int getDefaultConnectivity(Atom* a) {return getDefaultConnectivity(a->toString());}
	/**@brief Get the defaualt connectivity (not from forcefield data) for a given element*/
	int getDefaultConnectivity(const std::string& s)
	{
		for(const auto& p : defaultconnectivities)
		{
			if(get<0>(p)==s[0]) return get<1>(p);
		}
		cout << "Default Connectivity data not found for '"<<s[0]<<"'. Please enter: ";
		int v; cin >> v;
		defaultconnectivities.push_back(make_pair(s[0],v));
		return v;
	}
	/**@brief Get electronegativity data for a given element - experimental. If the data is not available, the user is prompted to enter it after which is is stored for that run*/
	double getElectronegativity(const std::string& s)
	{
		for(const auto& p : electronegativities)
		{
			if(get<0>(p)==s[0]) return get<1>(p);
		}
		cout << "Electronegativity data not found for '"<<s[0]<<"'. Please enter: ";
		double v; cin >> v;
		electronegativities.push_back(make_pair(s[0],v));
		return v;
	}
	inline double getScamGasteigerParameter(Atom* a) {return getElectronegativity(a->toString());}
	/**@brief Load additional atom-types - Useful when trying to add extra atom-types for specific usage (i.e. part of a molecule fragment).*/
	void loadAdditionalAtomTypes(const std::string& ff_file)
	{
		ifstream f;
    f.open(ff_file,ios::in);//opening file containing the data set
    int i=0;
    int v,h;
    double wt,c,s,e,r;
    while(true)
    {
      char *id=new char[5];
      f>>id>>v>>wt>>c>>s>>e>>r>>h;
      if(f.eof()) break;
      atom_types.push_back(new Atom(id,0,0,0,v,wt,c,s,e,r,h));
			delete[] id;
    }
	}
	/**@brief Load categories from an external file into this forcefield.
		 @details Categories can also be loaded when the ForceField object is created (See: ForceField(const std::string&,const std::string&))
		 @param[in] catf: Fileame to load category data from (See Atom-type group formats in documentation PDF)
		 @param[in] emptyCat: Empty existing category data? (Default: Yes)
	*/
	void loadCategories(const std::string& catf,bool emptyCat=true)
	{
		if(emptyCat) categories.empty();
		ifstream cf; cf.open(catf,ios::in);
		std::string cat,data;
		while(true)
		{
			cf >> cat >> data;
			if(cf.eof()) break;
			categories.push_back(Category(cat,data));
		}
		cout << "Categories loaded! "<<categories.size()<<" categor(y/ies) exist(s).\n";
	}
	/**@brief Load formal charges from a data file*/
	void loadChargeGroups(const std::string& chgfile)
	{
		ifstream cf; cf.open(chgfile,ios::in);
		std::string cat,data;
		while(true)
		{
			cf >> cat >> data;
			if(cf.eof()) break;
			Category tcat=Category(cat,data);
			for(Atom* a : atom_types) {if(tcat.contains(a)) a->setFormalCharge(std::stod(cat));}
		}
		cout << "Charge groups loaded!\n";
	}
	/**@brief Load solvation energy prediction weights from a file*/
	void loadSolvationWeights(const std::string& solvfile)
	{
		ifstream sf; sf.open(solvfile,ios::in);
		std::string cat,wt;
		while(true)
		{
			sf >> cat >> wt;
			if(sf.eof()) break;
			Category tcat=Category("temp",cat);
			for(Atom* a : atom_types)
			{
				if(tcat.contains(a)) a->setSolvationCorrection(std::stof(wt));
			}
		}
		cout << "Solvation weights loaded!\n";
	}
	/**@brief Load rules from a given file.
		 @details Rules define what atom-atom bonds are allowed, and how many. See the atom-type definition section of documentation PDF for more details on format.
		 @param[in] fl: Input filename
		 @param[in] emptyRules: Remove existing rule data? (Default: Yes)
	*/
	void loadRules(const std::string& fl,bool emptyRules=true) {rules=chemtools::getRuleset(fl,this); cout << "Loaded rules. Rules for "<<rules.size()<<" atom types found\n";}
	/**@brief Get the list of all the rules loaded into this forcefield
		 @details All the rules loaded into this ForceField are returned by atom-type.
		 @return A set of pairs of atom-type names and the set of rules imposed on it (std::vector<std::pair<std::string,std::vector<chemtools::Rule>>>)
	*/
	inline const auto& getRules() const {return rules;}
	/**@brief Checks if there are any rules defined on a given atom-type*/
	inline bool hasRule(const std::string& tp) const
	{
		for(const auto& d : rules)
		{
			//cout << get<0>(d)<<"\t"<<tp<<"\n";
			if(get<0>(d)==tp) return true;
		}
		return false;
	}
	inline bool hasRule(const Atom* at) const {return hasRule(at->toString());}
	/**@brief Checks if the given Category is present in the list loaded into this ForceField*/
	bool hasCategory(const std::string& name) const
	{
		for(const Category& c : categories) {if(c.getCategoryName()==name) return true;}
		return false;
	}
	/**@brief Gets the chosen Category by name if present in the list loaded into this ForceField. If there is no Category with this name, NoSuchCategoryException is thrown.*/
	const Category& getCategory(const std::string& name) const
	{
		for(const Category& c : categories) {if(c.getCategoryName()==name) return c;}
		throw NoSuchCategoryException();
	}
	/**@brief Uses the rules loaded into this ForceField to check if an atom (with the bonds supplied) satisfies all atom-type rules as defined (See Rule::satisfies(Atom*,Molecule*,const ForceField*))*/
	bool isSatisfied(Atom* a,const std::vector<Atom*>& bnd) const
	{
		if(a->mustCycl && !a->isCyclized()) return false;
		try
		{
			const auto& ruleset=findByKey(rules,a->toString());
			for(int i=0;i<ruleset.size();i++) {if(!ruleset[i].satisfies(a,bnd,this)) return false;}
			return true;
		}
		catch(ElementNotFoundException ex) {return true;}
	}
	/**@brief Uses the rules loaded into this ForceField to check if an atom (with the bonds supplied) satisfies all atom-type rules as defined given that an extra bond is formed to the new atom 'n' supplied (See Rule::satisfies(Atom*,Atom*,Molecule*))*/
	bool isSatisfied(Atom* a,Atom* n,std::vector<Atom*> bnd) const
	{
		//if(a->mustCycl && !a->isCyclized()) return false;
		try
		{
			const auto& ruleset=findByKey(rules,a->toString());
			for(int i=0;i<ruleset.size();i++) {if(!ruleset[i].satisfies(a,n,bnd)) return false;}
			return true;
		}
		catch(ElementNotFoundException ex) {return true;}
	}
	/**@brief Get a copy of an atom-type loaded into this ForceField by name
		 @details Returns a copy of the atom by name loaded with all force-field parameters. If no atom with such a name is loaded into the ForceField, a nullptr is returned
		 <br/>When using DeNovo with a well-parameterized Force-field, it is recommended that you use this method from a loaded "ForceField" object rather than directly using the Atom class' constructor.
	*/
  inline Atom* getAtom(const std::string& s) const
  {
    for(Atom* p : atom_types) {if(p->toString()==s) return new Atom(*p);}
    return nullptr;
  }
	/**@brief Same as getAtom(const std::string&) but uses an Atom object pointer instead of a string to get the atom-type name*/
  inline Atom* getAtom(Atom* a) const {if(a) return getAtom(a->toString()); else return getRandomNonHAtomByType();}
	/**@brief Randomly get a copy of any of the loaded atom-types*/
  inline Atom* getRandomAtomByType() const {return new Atom(*randomSelect(atom_types));}
	/**@brief Randomly get a copy of any of the loaded <b>non-hydrogen</b> atom-types*/
  inline Atom* getRandomNonHAtomByType() const
  {
	  Atom* t = randomSelect(atom_types); while(t->isHydrogen()) t=randomSelect(atom_types);
	  return new Atom(*t);
  }
	/**@brief Same as getRandomAtomByType() except that it also assignes the position to the atom at creation-time*/
  inline Atom* getRandomAtomByType(double x,double y,double z) const {return getAtom(getRandomAtomByType()->toString(),x,y,z);}
	/**@brief Same as getRandomNonHAtomByType() except that it also assignes the position to the atom at creation-time*/
  inline Atom* getRandomNonHAtomByType(double x,double y,double z) const {return getAtom(getRandomNonHAtomByType()->toString(),x,y,z);}
	/**@brief Get an atom of a specific atom-type at a specified position*/
  Atom* getAtom(const std::string& s,double x,double y,double z) const
  {
    int nu;
    Atom* natom=getAtom(s);
    //Atom* natom=new Atom(*tatm);
    natom->modify_x(x); natom->modify_y(y); natom->modify_z(z);
    return natom;
  }
	/**@brief Get the element name (symbol) corresponding to a given atom-type*/
	std::string getElementName(const std::string& el) const
	{
		if(hasCategory("halogen") && getCategory("halogen").contains(el))
		{
			switch(el[0])
			{
				case 'B': return "Br";
				case 'C': return "Cl";
			}
		}
		return std::string(1,el[0]);
	}
	/**@brief Get all atom-types corresponding to a given element*/
	std::vector<std::string> getAtomTypesByElement(const std::string& el)
	{
		std::vector<std::string> ret;
		const Category* halcat=nullptr;
		if(hasCategory("halogen")) halcat=&getCategory("halogen");
		if(halcat)
		{
			for(const std::string& ne : halcat->getGroup())
			{
				if(el.length()>1 || (el[0]!='C' && el[0]!='B')) {if(ne[0]==el[0] || ne[0]+32==el[0]) ret.push_back(ne);}
			}
			if(ret.size()) return ret;
		}
		for(const Atom* nep : atom_types)
		{
			if(halcat && halcat->contains(nep)) continue;
			if(nep->toString()[0]==el[0]) ret.push_back(nep->toString());
		}
		return ret;
	}
	/**@brief Get a random atom that can bond to this atom but bias the output
		 @details This function biases the output of the random result to return an atom-type of the same element as the biasing atom (actually it biases using the first character, so CL and BR cannot be bias targets)<br/>
		 If no such atom-type is available, then a random atom is returned. If multiple atom-types of the bias element are available, one is chosen randomly.
		 @param[in] a: The atom which has to form the bond
		 @param[in] bias: The atom-type (technically element) to which to bias the result
	*/
	Atom* selectRandomBondableAtom(Atom* a,Atom* bias=nullptr) const
	{
		std::vector<std::string> allatoms,biasedatoms;
		std::string temp=a->toString(),temp2;
		for(const BondData& bd : getAllowedBonds(a))
		{
			temp2=bd.getSecond(temp);
			allatoms.push_back(temp2);
			if(bias && temp2[0]==bias->toString()[0]) biasedatoms.push_back(temp2);
		}
		if(!allatoms.size()) throw NoBondsAvailableException();
		if(!biasedatoms.size()) return getAtom(randomSelect(allatoms));
		else return getAtom(randomSelect(biasedatoms));
	}
	/**@brief Get a random terminal atom that can bond to this atom but bias the output
		 @details This function biases the output of the random result to return an atom-type of the same element as the biasing atom (actually it biases using the first character, so CL and BR cannot be bias targets)<br/>If no such atom-type is available, then a random atom is returned. If multiple atom-types of the bias element are available, one is chosen randomly.<br/>
		 Note that this function specifically looks for a <b>terminal</b> atom. If no terminals are possible, a nullptr is returned.
		 @param[in] a: The atom which has to form the bond
		 @param[in] bias: The atom-type (technically element) to which to bias the result
	*/
	Atom* selectRandomBondableAtomTerminal(Atom* a,Atom* bias=nullptr) const
	{
		std::vector<Atom*> allatoms,biasedatoms;
		Atom* temp2=nullptr;
		for(const BondData& bd : getAllowedBonds(a))
		{
			temp2=bd.getSecond(a,*this);
			if(!temp2) continue;
			if(temp2->getStandardValency()==1)  allatoms.push_back(temp2);
			else {delete temp2; continue;}
			if(bias && temp2->toString()[0]==bias->toString()[0]) biasedatoms.push_back(temp2);
		}
		if(!allatoms.size()) return nullptr;
		Atom* reta;
		if(!biasedatoms.size()) reta=randomSelect(allatoms);
		else reta=randomSelect(biasedatoms);
		for(Atom*& ap : allatoms) {if(ap!=reta) delete ap;}
		return reta;
	}
	/**@brief Situationally list atoms that can bond to the selected atom
		 @details Get a list  of the loaded atom-types such that they can bond to the source atom (see params).<br/>
		 This function ensures that after adding this bond, the no atom-type definitions are violated. No checks are made to ensure that all the rules can still be satisfied after this stage. It is assumed that the rules are sane enough to allow the atom-type to satisfy them
	*/
	inline std::vector<Atom*> listAtomsByRule(Atom* s,Molecule* m,Atom* temp=nullptr) const;
	/**@brief Same as listAtomsByRule(Atom*,Molecule*,Atom*), but directly uses the list of bonded atoms instead of collecting them from the Molecule object*/
	std::vector<Atom*> listAtomsByRule(Atom* s,const std::vector<Atom*>& bnd,Atom* temp=nullptr) const
	{
		if(!rules.size()) cout << "WARN: Ruleset not loaded?\n";
		std::vector<chemtools::Rule> ruleset;
		try
		{
			ruleset=findByKey(rules,s->toString());
			int ind=-1;
			#ifdef SCHUFFLE_RULES
			std::vector<chemtools::Rule> allrules;
			#endif
			for(int i=0;i<ruleset.size();i++)
			{
				if(ruleset[i].satisfies(s,bnd)) continue;
				ind=i;
				#ifdef SCHUFFLE_RULES
				allrules.push_back(ruleset[i]);
				#else
				break;
				#endif
			}
			if(ind!=-1)
			{
				#ifdef SCHUFFLE_RULES
				return std::vector<Atom*>(1,randomSelect(allrules).suggest(s,this));
				#else
				return std::vector<Atom*>(1,ruleset[ind].suggest(s,this));
				#endif
			}
		}
		catch(ElementNotFoundException ex) {}
		std::vector<Atom*> allowedatoms;
		//Atom* ret=randomSelect(getAllowedBonds(s)).getSecond(s,*this);
		bool sat=true;
		for(const BondData&  b: getAllowedBonds(s))
		{
			Atom* ret=b.getSecond(s,*this);
			if(!ret) continue;
			sat=true;
			for(chemtools::Rule& r : ruleset) {if(!r.satisfies(s,ret,bnd)) {sat=false; break;}}
			if(sat) allowedatoms.push_back(ret);
			else delete ret;
		}
		if(!allowedatoms.size()) throw RuleConstraintsPreventBondsException();
		else
		{
			if(temp)
			{
				char chk=temp->toString()[0];
				std::vector<Atom*> f2;
				std::vector<Atom*> refats;
				for(Atom* a : allowedatoms) {if(a->toString()[0]==chk) f2.push_back(a); else refats.push_back(a);}
				if(f2.size()) {allowedatoms=f2; for(Atom* tda : refats) delete tda;}
			}
			return allowedatoms;
		}
	}
	/*@brief Select a random atom from the list of possible bondable atoms (See listAtomsByRule(Atom*,Molecule*,Atom*))*/
	inline Atom* selectAtomByRule(Atom* s,Molecule* m,Atom* temp=nullptr) const;
	/**@brief Same as selectAtomByRule(Atom*,Molecule*,Atom*), but directly uses the list of bonded atoms instead of collecting them from the Molecule object*/
	Atom* selectAtomByRule(Atom* s,const std::vector<Atom*>& bnd,Atom* temp=nullptr) const
	{
		auto allowedatoms=listAtomsByRule(s,bnd,temp);
		Atom* r=randomSelect(allowedatoms);
		for(Atom*& ap : allowedatoms) {if(ap!=r) delete ap;}
		return r;
	}
	/**@brief List all possible a <b>fragments</b> that can bond to a given atom (following bonding rules).
		 @details The output also contains the target atom and its "support" (two atoms connected to it)
	*/
	inline std::vector<std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>>> listFragmentsByRule(Atom* src,Molecule* m,const FragmentSet& frags,Atom* exa=nullptr,int devid=0) const;
	/**@brief Same as listFragmentsByRule(Atom*,Molecule*,const FragmentSet&,Atom*), but directly uses the list of bonded atoms instead of collecting them from the Molecule object*/
	std::vector<std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>>> listFragmentsByRule(Atom* src,const std::vector<Atom*>& bnd,const FragmentSet& frags,Atom* exa=nullptr,int devid=0) const;
	/*@brief Select a random <i>allowed</i> fragment from the list (See listFragmentsByRule(Atom*,Molecule*,const FragmentSet&,Atom*))*/
	inline std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> selectFragmentByRule(Atom* s,Molecule* m,const FragmentSet& fs,Atom* temp=nullptr) const;
	/**@brief Same as selectFragmentByRule(Atom*,Molecule*,const FragmentSet&,Atom*), but directly uses the list of bonded atoms instead of collecting them from the Molecule object*/
	inline std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> selectFragmentByRule(Atom* s,const std::vector<Atom*>& bnd,const FragmentSet& fs,Atom* temp=nullptr) const
	{
		auto opts=listFragmentsByRule(s,bnd,fs,temp);
		if(opts.size()) return randomSelect(opts);
		else throw RuleConstraintsPreventBondsException();
	}
	/**@brief Checks whether accepting a particular atom to bond to a source atom is allowed accornding to atom-type definitions.
		 @details This function is similar to selectAtomByRule(), except instead of providing the atom to be chosen, it checks if the chosen atom satisfies the required definitions.<br/>
		 It is used when checking for cyclization
		 @param[in] src: The source atom on which the check must be performed
		 @param[in] bnd: The list of all atoms already bonded to src
		 @param[in] n: The new atom being trialled for bonding to src
	*/
	bool acceptNewBond(Atom* src,std::vector<Atom*> bnd,Atom* n) const
	{
		if(!rules.size()) cout << "WARN: Ruleset not loaded?\n";
		std::vector<chemtools::Rule> ruleset;
		try
		{
			ruleset=findByKey(rules,src->toString());
			for(int i=0;i<ruleset.size();i++)
			{
				if(ruleset[i].satisfies(src,n,bnd)) continue;
				else return false;
			}
			return true;
		}
		catch(ElementNotFoundException ex) {return true;}
	}

	/**@brief See contains(std::string)*/
  bool contains(Atom* a) const
  {
    for(Atom* p : atom_types)
    {
      if(p==a)return true;
    }
    return false;
  }
	/**@brief Checkes if this ForceField has an atom-type with the given name*/
  bool contains(std::string str) const
  {
    for(Atom* p : atom_types)
    {
      if(p->toString()==str)return true;
    }
    return false;
  }
	/**@brief Write valency, σ, and ε of all the atom-types to stdout*/
  void dump()
  {
    for(Atom* a : atom_types)
      cout << a->toString()<<","<<a->seek_valency()<<","<<a->seek_sigma()<<","<<a->seek_epsilon()<<"\n";
  }
	/**@brief Get a list of all the bonds (as BondData objects) that this atom is allowed to make.
		 @details Note that this is different from selectAtomByRule() which considers the atoms bonded to the source atom.<br/>This function only provides a list of all possible bonds that this atom-type can have
	*/
  inline const std::vector<BondData>& getAllowedBonds(Atom* a) const {return getAllowedBonds(a->toString());}
	/**@brief Same as getAllowedBonds(Atom*), but provides the atom-type by name (as std::string)*/
  const std::vector<BondData>& getAllowedBonds(const std::string& temp) const
  {
    for(int i=0;i<atombondmap.size();i++) {if(get<0>(atombondmap[i])==temp) return get<1>(atombondmap[i]);}
    throw NoBondsAvailableException();
  }
	/**@brief Expand dihedral multiplicity to explicitly write down the allowed dihedrals. See algorithm PDF (the dihedrals section, and attached links) for more details.*/
	static std::vector<double> expandMultiplicity(double v,int m) //v is in radians
	{
		double EX=(2*PI)/m;
		std::vector<double> ret;
		for(int i=0;i<m;i++) ret.push_back(v+i*EX);
		return ret;
	}

	/**@brief Load bond-length data from a given file
		 @param[in] s: Filename to load from.
		 @param[in] erase: Erase any pre-existing data (Default: Yes)
	*/
  void loadLengths(const std::string& s,bool erase=true);
	/**@brief Load bond-angle data from a given file
		 @param[in] s: Filename to load from.
		 @param[in] erase: Erase any pre-existing data (Default: Yes)
	*/
  void loadAngles(const std::string& s,bool erase=true);
	/**@brief Load dihedral data from a given file
		 @param[in] s: Filename to load from.
		 @param[in] erase: Erase any pre-existing data (Default: Yes)
	*/
  void loadDihedrals(const std::string& s,bool erase=true);
	/**@brief Load improper dihedral data from a given file (not used in later versions of DeNovo)
		 @details DeNovo Versions after SIGMA do not deal with improper dihedrals except for fallback cases or for programming explicit requirements.<br/>
		 The data for those parameters has not been updated since, and many new atom types (and parameters for existing atom-types) have been added for all other bonding data, so be sure to check the data for completeness and validity before using this feature.<br/>
		 <b>Warning: </b> For all practical purposes and general use of the DeNovo, this feature is considered <i>deprecated</i>.
	*/
  void loadImpropers(const std::string& s);
	/**@brief A wrapper to load all parameters (bond-length,bond-angle,dihedral,and improper-dihedral). See documentation (PDF) for default filenames
		 @details This function is a wrapper to load all required parameters at once using default (or, based on programmer's choice, other specified) input files for the data.<br/>
		 To load any of bond length, angle, dihedral, or improper parameters separately, look at the respective functions (ForceField::loadLengths, ForceField::loadAngles, ForceField::loadDihedrals, ForceField::loadImpropers).
		 @param[in] lengthfn: Filename to load <u>length</u> data (1st parameter)
		 @param[in] anglefn: Filename to load <u>angle</u> data (2nd parameter)
		 @param[in] dihedralfn: Filename to load <u><b>proper</b> dihedral</u> data (3rd parameter)
		 @param[in] improperfn: Filename to load <u><b>improper</b> dihedral</u> data (4th parameter)
		 @param[in] erase: Delete any pre-existing data? (Default: Yes)
	*/
  void loadBondParameters(std::string lengthfn=LengthDataFile,std::string anglefn=AngleDataFile,std::string dihedralfn=DihedralDataFile,std::string improperfn=ImproperDataFile,bool erase=true)
  {
    //cout << lengthfn << "; Now working\n";
    loadLengths(lengthfn,erase);
    loadAngles(anglefn,erase);
    loadDihedrals(dihedralfn,erase);
    loadImpropers(improperfn);
    cout << "Done\n";
  }
	/**@brief Load one-four interaction parameters from a given file. Later versions of DeNovo directly use the parameters from the atom-type data*/
  void loadOneFourParameters(const std::string&,const std::string&);
	/**@brief Load L-J parameters from a given file. Later versions of DeNovo use the Atom object for parameters directly*/
  void loadOtherInteractions(const std::string&);
	/**@brief Load 1-4 non-bonding potential parameters, as well as L-J potential parameters*/
  void loadNonbondedInteractions(const std::string& onefourfile=OneFourFile,const std::string& otherfile=OtherInterFile)
  {
    loadOneFourParameters(onefourfile,otherfile);
    loadOtherInteractions(otherfile);
  }

	/**@brief Get the bond length for the bond between two atom types.
		 @details Searches the bond data loaded into this force-field to find bond parameters for a given pair of atom types.<br/>
		 If it is not found (these two atom-types do not have a bond parameter), it will print an error message (unless silent mode is turned on), and throw a DataNotAvailableException
	*/
	inline double getBondLength(Atom* a1,Atom* a2,bool silent=false) const {return get<0>(getLengthParameters(a1->toString(),a2->toString(),silent));}
	/**@brief Alternate method for getLengthParameters(Atom*,Atom*,bool) using strings instead of Atom pointers*/
  std::pair<double,double> getLengthParameters(const std::string& a1,const string& a2,bool silent=false) const {return getLengthParameters(make_pair(a1,a2),silent);}
	/**@brief Alternate method for getLengthParameters(Atom*,Atom*,bool) using string lists instead of Atom pointers*/
  inline std::pair<double,double> getLengthParameters(const std::vector<std::string>& sv,bool silent=false) const {return getLengthParameters(sv[0],sv[1],silent);}
	/**@brief Alternate method for getLengthParameters(Atom*,Atom*,bool) using string pairs instead of Atom pointers*/
  inline std::pair<double,double> getLengthParameters(const std::pair<std::string,std::string>& atp,bool silent=false) const
  {
    for(auto& a : length_parameters)
    {
      if(matches<std::string>(a.getAtoms(),atp))
        return (a.getParameters());
    }
		if(!silent)
		{
	    cout << "Length data not available\n";
			cout << get<0>(atp)<<","<<get<1>(atp) << "\n";
		}
    //return make_pair(0,0); //Not found
    throw DataNotAvailableException();
  }
	/**@brief Get the bond angle for the angle between three atom types. -failsafe method
		 @details Searches the angle data loaded into this force-field to find angle parameters for a given triplet of atom types. It also checks if any templates match the request (see input files section of documentation PDF)<br/>
		 If it is not found, it will print a warning message (stating that data could not be found), and will resort to a default value (instead of throwing an exception - which is why this method has "failsafe" in the name).
	*/
	inline const std::vector<double>& getAngleParamtersFailsafe(Atom* a1,Atom* a2,Atom* a3) const {return getAngleParamtersFailsafe(a1->toString(),a2->toString(),a3->toString());}
	/**@brief Alternate method for getAngleParamtersFailsafe(Atom*,Atom*,Atom*) using strings instead of Atom pointers*/
  const std::vector<double>& getAngleParamtersFailsafe(const std::string& s1,const std::string& s2, const std::string& s3) const
  {
    try{return getAngleParamters(s1,s2,s3);}
    catch(DataNotAvailableException ex)
    {
      try{return getAngleParamters("X",s2,"X");}
      catch(DataNotAvailableException ex)
      {
        std::string t1(1,s1[0]),t2(1,s2[0]),t3(1,s3[0]);
        try {return getAngleParamters(t1,t2,t3);}
        catch(DataNotAvailableException ex)
        {
          try {return getAngleParamters("X",t2,"X");} catch(DataNotAvailableException ex) {cout << "Angle data not available ("<<s1<<","<<s2<<","<<s3<<"). Resorting to default\n"; return getAngleParamters("DUM","DUM","DUM");}
        }
      }
    }
  }
	/**@brief Alternate method for getAngleParamtersFailsafe(Atom*,Atom*,Atom*) using a string list instead of Atom pointers*/
  inline const std::vector<double>& getAngleParamtersFailsafe(const std::vector<std::string>& ats) const {return getAngleParamtersFailsafe(ats[0],ats[1],ats[2]);}
	/**@brief Get the bond angle for the angle between three atom types.
		 @details Searches the angle data for a <b>direct match</b> only (no template matching). If angle data is not found, DataNotAvailableException is thrown.
	*/
  const std::vector<double>& getAngleParamters(const std::string& a1,const std::string& a2, const std::string& a3) const
  {
    std::vector<std::string> atl; atl.push_back(a1);  atl.push_back(a2); atl.push_back(a3);
    return getAngleParamters(atl);
  }
	/**Alternate method for getAngleParamters(const std::string&,const std::string&,const std::string&) with all strings in a single list*/
  const std::vector<double>& getAngleParamters(const std::vector<std::string>& ats) const
  {
    std::vector<double> par; //th0,cth,ub0,cub
    int lev=4;
    std::vector<std::string> fb1,fb2,fb3;
    for(const auto& a : angle_parameters) {if(matches<std::string>((get<0>(a)),ats)) {return (get<1>(a));}}
    throw DataNotAvailableException();
    /*if(par.size()) return par;
    else {if(fallback) return getAngleParamters("C","C","C",false); else throw DataNotAvailableException();} //Get defaults and save later*/
  }

	/**@brief Get the dihedral angle for a set of four atoms (in order of the chain). -failsafe method
		 @details Searches the dihedral data loaded into this force-field to find dihedral parameters for a given set of atom types. It also checks if any templates match the request (see input files section of documentation PDF)<br/>
		 If it is not found, it will print a warning message (stating that data could not be found), and will resort to a default value (instead of throwing an exception - which is why this method has "failsafe" in the name).
	*/
  inline const std::vector<double>& getDihedralParametersFailsafe(Atom* a1,Atom* a2,Atom* a3,Atom* a4) const {return getDihedralParametersFailsafe(a1->toString(),a2->toString(),a3->toString(),a4->toString());}
	/**@brief Alternate method for getDihedralParametersFailsafe(Atom*,Atom*,Atom*,Atom*) using a strings instead of Atom pointers*/
  const std::vector<double>& getDihedralParametersFailsafe(const std::string& s1,const std::string& s2,const std::string& s3,const std::string& s4) const
  {
    std::vector<double> par; //th0,cth,ub0,cub
    try{return getDihedralParameters(s1,s2,s3,s4);}
    catch(DataNotAvailableException ex)
    {
			try{return getDihedralParameters("X",s2,s3,"X");}
      catch(DataNotAvailableException ex)
      {
				std::string t1(1,s1[0]),t2(1,s2[0]),t3(1,s3[0]),t4(1,s4[0]);
	      try {return getDihedralParameters(t1,t2,t3,t4);}
        catch(DataNotAvailableException ex)
        {
          try {return getDihedralParameters("X",t2,t3,"X");} catch(DataNotAvailableException ex) {cout << "Dihedral data not available("<<s1<<","<<s2<<","<<s3<<","<<s4<<"). Resorting to default\n"; return getDihedralParameters("C","C","C","C");}
        }
      }
    }
  }
	/**@brief Get the dihedral angle data for a set of four atom-types (in order of the connection chain)
		 @details Searches the dihedral data for a <b>direct match</b> only (no template matching). If dihedral data is not found, DataNotAvailableException is thrown.
	*/
  const std::vector<double>& getDihedralParameters(const std::string& a1,const std::string& a2,const std::string& a3,const std::string& a4) const //,bool fallback=true) const
  {
    std::vector<std::string> atl; atl.push_back(a1);  atl.push_back(a2); atl.push_back(a3); atl.push_back(a4);
    return getDihedralParameters(atl); //,fallback);
    //return getDihedralParameters(atl);
  }
	/**Alternate method for getDihedralParameters(const std::string&,const std::string&,const std::string&,const std::string&) with all strings in a single list*/
  const std::vector<double>& getDihedralParameters(const std::vector<std::string>& ats) const //,bool fallback=true) const
  {
    for(const auto& a : dihedral_parameters)
    {
      if(matches<std::string>((get<0>(a)),ats))
        return (get<1>(a));
    }
    //cout << "Dihedral parameters are not available for dihedral with:\t";
    //for(const std::string& s : ats) {cout <<s<<" ";}
    throw DataNotAvailableException();
  }

	/**@brief Exactly matches the atom-type names to the improper dihedral data and retrieves the angle. If not available DataNotAvailableException is thrown*/
  const std::pair<double,double>& getImproperParameters(const std::string& a1,const std::string& a2,const std::string& a3,const std::string& a4) const
  {
    std::vector<std::string> atl; atl.push_back(a1);  atl.push_back(a2); atl.push_back(a3); atl.push_back(a4);
    return getImproperParameters(atl);
  }
	/**@brief Alternate method for getImproperParameters(const std::string&,const std::string&,const std::string&,const std::string&) with all strings in a single list*/
  const std::pair<double,double>& getImproperParameters(const std::vector<std::string>& ats) const
  {
    for(auto& a : improper_parameters)
    {
      if(matches<std::string>((get<0>(a)),ats))
        return (get<1>(a));
    }
    //cout << "Improper parameters not available\n";
    throw DataNotAvailableException();
  }

	const std::pair<double,double> getImproperParametersFailsafe(const std::vector<std::string>& ats,int nb,int mv) const
	{
		try{return getImproperParameters(ats);}
		catch(DataNotAvailableException ex)
		{
			if(nb<=1 || mv<=2) return make_pair(0,0);
			if(nb==2) return make_pair(0,400);
			else return make_pair(0.59,60);
		}
	}
	inline const std::pair<double,double> getImproperParametersFailsafe(const std::string& a1,const std::string& a2,const std::string& a3,const std::string& a4,int nb,int mv) const
	{
		std::vector<std::string> ats(1,a1); ats.push_back(a2); ats.push_back(a3); ats.push_back(a4);
		return getImproperParametersFailsafe(ats,nb,mv);
	}
	inline const std::pair<double,double> getImproperParametersFailsafe(Atom* a1,Atom* a2,Atom* a3,Atom* a4) const
	{
		std::vector<std::string> ats(1,a1->toString()); ats.push_back(a2->toString()); ats.push_back(a3->toString()); ats.push_back(a4->toString());
		return getImproperParametersFailsafe(ats,a1->getHybridization(),a1->getStandardValency());
	}

  inline const std::pair<double,double>& getOneFourParameters(const std::string& a1,const std::string& a2) const {return getOneFourParameters(make_pair(a1,a2));}
  const std::pair<double,double>& getOneFourParameters(const std::pair<std::string,std::string>& v) const
  {
    for(auto& a : one_four_interactions)
    {
      if(matches((get<0>(a)),v))
        return (get<1>(a));
    }
    //cout << "One-Four parameters not available\n";
    throw DataNotAvailableException();
    //cout << "Not found(One-Four Parameters)\n";
  }
	/**@brief Get the σ and ε values for a given atom-type*/
  inline const std::vector<double>& getOtherParameters(const std::string& s) const {return getCorresponding(other_interactions,s);}

  //Calculating the energies
	/**@brief Calculate the angle-strain potential energy (see documentation in PDF for the exact function)
		 @details Calculates the net energy from angle parameters (including Urey-Bradley contribution - See PDF)<br/>This method can also be used as a boolean tolerance (with the help of the "rigid" flag)
		 @param[in] rigid: A boolean flag. If set to true, instead of calculating the energy, the function checks if the angle value falls within a tolerance limit (acut) of the expected value. If yes, the energy is taken to be zero, otherwise it's ∞
		 @param[in] acut: The angle tolerance value. The default tolerance is THETACUT
	*/
  double calculateAngleEnergy(Atom* a0,Atom* a1,Atom* a2,bool rigid=false,double acut=THETACUT) const
  {
		double cv=(((a0->seek_x())-(a1->seek_x()))*((a2->seek_x())-(a1->seek_x()))+((a0->seek_y())-(a1->seek_y()))*((a2->seek_y())-(a1->seek_y()))+((a0->seek_z())-(a1->seek_z()))*((a2->seek_z())-(a1->seek_z())))/(sqrt(pow(((a0->seek_x())-(a1->seek_x())),2)+pow(((a0->seek_y())-(a1->seek_y())),2)+pow(((a0->seek_z())-(a1->seek_z())),2))*sqrt(pow(((a2->seek_x())-(a1->seek_x())),2)+pow(((a2->seek_y())-(a1->seek_y())),2)+pow(((a2->seek_z())-(a1->seek_z())),2)));
		double th=0;
		if(cv>=1) th=0;
		else if(cv<=-1) th=PI;
    else th=acos(cv);
    return calculateAngleEnergy(a0->toString(),a1->toString(),a2->toString(),th,a2->distanceFrom(a0),rigid,acut);
  }
	/**@brief Alternate method for calculateAngleEnergy(Atom*,Atom*,Atom*,bool,double) using strings instead of Atom pointers*/
  double calculateAngleEnergy(const std::string& s1,const std::string& s2,const std::string& s3,double th,double ub,bool rigid=false,double acut=THETACUT) const
  {
    std::vector<double> par=getAngleParamtersFailsafe(s1,s2,s3); //th0,cth,ub0,cub
    if(!par.size()) throw DataNotAvailableException(); //{return 0;} //Data not available
    if(rigid) {if(abs(par[0]-th)>acut) {return 1e20;} else {return 0;}}
    else
    {
      //par[0]*=0.01745329252; (Corrected while loading itself)
			return angleEnergyFunction(par,th,ub);
      //return angle_energy+urey_bradley_energy;
    }
  }

	/**@brief Calculate the dihedral potential energy (see documentation in PDF for the exact function)
		 @details Calculates the net energy from dihedral parameters.<br/>This method can also be used as a boolean tolerance (with the help of the "rigid" flag)
		 @param[in] rigid: A boolean flag. If set to true, instead of calculating the energy, the function checks if the dihedral angle value falls within a tolerance limit (acut) of the expected value. If yes, the energy is taken to be zero, otherwise it's ∞
		 @param[in] acut: The dihedral tolerance value. The default tolerance is DIHEDTHETACUT
	*/
  double calculateDihedralEnergy(Atom* a1,Atom* a2,Atom* a3,Atom* a4,bool rigid=false,double acut=DIHEDTHETACUT) const
  {
    double phi=quickgeom::dihedralBetween(a1->bondVectorTo(a2),a2->bondVectorTo(a3),a3->bondVectorTo(a4));
    return calculateDihedralEnergy(a1->toString(),a2->toString(),a3->toString(),a4->toString(),phi,rigid,acut);
  }
	/**@brief Alternate method for calculateDihedralEnergy(Atom*,Atom*,Atom*,Atom*,bool,double) using strings instead of Atom pointers*/
  double calculateDihedralEnergy(const std::string& s1,const std::string& s2,const std::string& s3,const std::string& s4,double ph,bool rigid=false,double acut=DIHEDTHETACUT) const
  {
    const std::vector<double>& par=getDihedralParametersFailsafe(s1,s2,s3,s4); // phi0, K, n
    if(!par.size()) throw DataNotAvailableException(); //{return 0;} //Data not available
    if(rigid) {if(abs(par[0]-ph)>acut) {return 1e15;} else {return 0;}}
    else return dihedralEnergyFunction(par,ph);
  }
	/**@brief Calculate inproper dihedral energy. Not used in latest versions of DeNovo. Refer to links provided in documentation PDF. Throws DataNotAvailableException if data is not available*/
  inline double calculateImproperEnergy(Atom* a0,Atom* a1,Atom* a2,Atom* a3) const
  {
    double *m=cross_product(a1->seek_x()-a0->seek_x(),a1->seek_y()-a0->seek_y(),a1->seek_z()-a0->seek_z(),a2->seek_x()-a0->seek_x(),a2->seek_y()-a0->seek_y(),a2->seek_z()-a0->seek_z()),*n=cross_product(a2->seek_x()-a0->seek_x(),a2->seek_y()-a0->seek_y(),a2->seek_z()-a0->seek_z(),a3->seek_x()-a0->seek_x(),a3->seek_y()-a0->seek_y(),a3->seek_z()-a0->seek_z());
		double cv=dot_product(m[0],m[1],m[2],n[0],n[1],n[2])/(sqrt(pow(m[0],2)+pow(m[1],2)+pow(m[2],2))*sqrt(pow(n[0],2)+pow(n[1],2)+pow(n[2],2)));
		double q=0;
		if(cv>=1) q=0;
		else if(cv<=-1) q=PI;
    else q=acos(cv);
    delete[] m; delete[] n;
    return calculateImproperEnergy(a0->toString(),a1->toString(),a2->toString(),a3->toString(),q);
  }
	/**@brief Alternate method for calculateImproperEnergy(Atom*,Atom*,Atom*,Atom*) using strings instead of Atom pointers.*/
  double calculateImproperEnergy(const std::string& s1,const std::string& s2,const std::string& s3,const std::string& s4,const double& q) const
  {
    try{std::pair<double,double> par=getImproperParameters(s1,s2,s3,s4); //q0,cq
    double q0=get<0>(par),cq=get<1>(par);
    //q0*=0.01745329252;
    return ((q0-q)*(q0-q)*cq)/2;}
    catch(DataNotAvailableException ex) {return 0;}
  }
};


void collateDataFiles(const std::string& atomtype_file="atomtypes.atp",const std::string& ljc_file="ljcparams.dnvin",const std::string& outfile="final_ff_parameters.ffin")
{
  ifstream f;
  std::vector<std::pair<std::string,int>> valences;
  std::vector<std::pair<std::string,double>> charges;
  std::vector<std::pair<std::string,double>> weights;
  std::vector<std::pair<std::string,double>> radii;
  std::vector<std::pair<std::string,std::pair<double,double>>> sigmaepsilons;
  std::vector<std::string> atomtypes;

  std::string attp;
  double d1,d2,d3;
  int i1;

  //Loading weights and valences from atomtypes.atp (of ff)
  cout << "Loading Mass and Valences\n";
  f.open(atomtype_file,ios::in);//opening file containing the data set
  if(!f) cout << "File not opened\n";
  while(true)
  {
    f >> attp >> d1>>i1;
    if(f.eof()) break;
    atomtypes.push_back(attp);
    valences.push_back(make_pair(attp,i1));
    weights.push_back(make_pair(attp,d1));
  }
  f.close();
  cout <<"Done\n";

  //Loading the LJ Parameters from (dnv) file
  cout << "Loading LJ Parameters\n";
  f.open(ljc_file,ios::in);//opening file containing the data set
  while(true)
  {
    f >> attp >> d3>>d1>>d2;
    cout << attp <<" "<< d3 << " "<<d1 << " "<< d2<<"\n";
    if(f.eof()) break;
    //atomtypes.push_back(attp);
    sigmaepsilons.push_back(make_pair(attp,make_pair(d1,d2)));
    charges.push_back(make_pair(attp,d3));
  }
  f.close();
  cout <<"Done\n";

  //Loading the Charge parameters

  ofstream wo;
  wo.open(outfile,ios::out);
  //f>>id>>v>>wt>>c>>s>>e>>r;
  std::pair<double,double> ljp;
  for(std::string& s : atomtypes)
  {
    ljp=getCorresponding(sigmaepsilons,s);
    if(!(get<0>(ljp))) continue;
    wo << s <<" "<<getCorresponding(valences,s)<<" "<<getCorresponding(weights,s)<<" "<<getCorresponding(charges,s)<<" "<<(get<0>(ljp))<<" "<<(get<1>(ljp))<<" "<<(get<0>(ljp))/2<<"\n"; // !!!NOTE: sigma/2 is the vanderwaal_radius
  }
  wo.close();
  cout << "Written to: "<<outfile<<"\n";
}

//Completing the Bond class
inline Atom* BondData::getSecond(Atom* a,const ForceField& ff) const { return ff.getAtom(getSecond(a->toString()));}

//Completing Category class
Category::Category(const std::string& s1,const std::string& rl,const ForceField* ff)
{
	targa=s1;
	reps=rl;
	std::vector<std::string> tempgroup=stringfx::split(rl,'|',true);
	myGroup=std::vector<std::string>();
	if(ff)
	{
		for(std::string s : tempgroup)
		{
			for(Atom* a : ff->getAtomTypes()) {if(stringfx::matches(s,a->toString())) myGroup.push_back(a->toString());}
		}
	}
	else myGroup=tempgroup;
}
static std::ostream& operator<<(ostream& os,const Residue& res)
{
	os << res.getName()<<": ";
	for(int i=0;i<res.getAtomCount();i++) os << res.getAlias(i)<<"("<<res.getActual(i)<<")-";
	return os;
}
/**@brief Load residues from an RTP file
	 @details Standard RTP files can be used by DeNovo for loading residues. Each residue has unique atom names for each atom which is then mapped to a forcefield atom type. (Open protein.rtp attached with DeNovo to see the format)<br/>
	 These files can be found in GROMACS data files, or (with enough effort) written manually. It's preferred to use GROMACS rtp files so that GROMACS can be used to convert/rewrite PDB files (and hence the atom-type names and unique names of all the atoms match)<br/>
	 The protein.rtp and dna.rtp files are taken from GROMACS data files for the CHARMM-27 forcefield based on which the DeNovo	generation is based as well.
*/
void ForceField::loadResidues(const std::string& resf)
{
	ifstream fs; fs.open(resf,ios::in);
	std::string curline,section,resname;
	std::vector<std::string> lineset;
	residues=std::vector<Residue>();
	bool appmode=false;
	while(true)
	{
		getline(fs,curline);
		if(fs.eof()) break;
		int bI=stringfx::indexOf('[',curline);
		if(bI!=-1)
		{
			int ebI=stringfx::indexOf(']',curline);
			section=stringfx::trim(curline.substr(bI+1,ebI-bI-1));
			if(section=="atoms") {appmode=true; continue;}
			else
			{
				appmode=false;
				if(lineset.size())
				{
					residues.push_back(Residue(resname,lineset));
					lineset=std::vector<std::string>();
				}
				resname=section;
				continue;
			}
		}
		if(appmode) lineset.push_back(curline);
	}
	cout << residues.size()<<" residues loaded.\n";
}
/**\class IncompleteForceField forcefield.hpp "forcefield/forcefield.hpp"
	 A force-field class to load force-field with missing valency and hybridization data. It used to load atoms into macromolecules (where bond calculation is unnecessary)
*/
class IncompleteForceField : public ForceField
{
public:
	/**@brief Load an incomplete force-field from a given file - The default IncompleteForceField class constructor*/
	IncompleteForceField(const std::string& ff_file) : ForceField()
  {
    atom_types=std::vector<Atom*>();
    ifstream f;
    f.open(ff_file,ios::in);//opening file containing the data set
    cout << "Loading (Incomplete) Forcefield from: "<<ff_file<<"\n";
    int i=0;
    int v,h;
    double wt,c,s,e,r;
    while(true)
    {
      char *id=new char[5];
      f>>id>>wt>>c>>s>>e;
			r=s/2;
      if(f.eof()) break;
      atom_types.push_back(new Atom(id,0,0,0,4,wt,c,s,e,r));
			delete[] id;
    }
    cout <<"Incomplete ForceField loaded! "<<atom_types.size()<<" atom types available and have some (missing valency and hybridization) necessary parameters.\n";
  }
};

//Completing the chemtools namespace
chemtools::Rule::Rule(std::string src,std::string rule,const ForceField* ff)
{
	centre=src;
	std::vector<std::string> pcs=stringfx::split(rule,' ',true);
	group=stringfx::split(pcs[0],'|');
	try
	{
		if(ff)
		{
			std::vector<std::string> fgroup;
			for(std::string& s : group) for(const BondData& b : ff->getAllowedBonds(centre)) if(stringfx::matches(s,b.getSecond(centre))) fgroup.push_back(b.getSecond(centre));
			if(fgroup.size()) group=fgroup;
		}
	}
	catch(NoBondsAvailableException ex) {cout << "No bonds for atom type: "<<centre<<"\n";}
	//cout << pcs << "\n";
	min=std::stoi(pcs[1]);
	if(pcs.size()>2 && pcs[2][0]!='-') max=std::stoi(pcs[2]);
	if(pcs.size()>3 && pcs[3][0]!='-') {angp=true; angv=stod(pcs[3]);}
	if(pcs.size()>4 && pcs[4][0]!='-')
	{
		dihp=true;
		std::vector<std::string> pcsd=stringfx::split(pcs[4],',');
		dihv=stod(pcsd[0]);
		dihm=stod(pcsd[1]);
	}
}
std::vector<chemtools::Rule> chemtools::buildRules(const std::string& def,const ForceField* ff)
{
	int i=stringfx::indexOf(' ',def);
	std::string c=def.substr(0,i);
	std::string r=def.substr(i+1,def.length()-(i+1));
	std::vector<std::string> rstr=stringfx::split(r,';');
	std::vector<chemtools::Rule> rules;
	for(std::string& s : rstr)
	{
		if(ff) for(const auto& p : ff->getCategories()) stringfx::replace(s,"["+p.getCategoryName()+"]",p.toString());
		rules.push_back(chemtools::Rule(c,s,ff));
	}
	return rules;
}
std::vector<std::pair<std::string,std::vector<chemtools::Rule>>> chemtools::getRuleset(const std::string& file,const ForceField* ff)
{
	std::vector<std::pair<std::string,std::vector<chemtools::Rule>>> ret;
	std::ifstream inf; inf.open(file,ios::in);
	if(!inf) throw FileNotFoundException();
	std::string str;
	int constcount=0;
	while(true)
	{
		getline(inf,str);
		if(inf.eof()) break;
		std::vector<chemtools::Rule> rules=buildRules(str,ff);
		if(!rules.size()) continue;
		if(ff)
		{
			//cout << rules[0].getCentre()<<" "<< rules.size() <<"\n";
			constcount=0;
			for(const Rule& r : rules) constcount+=r.getMin();
			#ifndef QUIET
			cout << rules[0].getCentre()<<" "<< constcount <<"\n";
			#endif
			//if(constcount>ff->getAtom(rules[0].getCentre())->seek_valency()) throw RulesExceedValenciesException();
		}
		ret.push_back(make_pair(rules[0].getCentre(),rules));
	}
	return ret;
}
inline std::vector<double> chemtools::Rule::predictAngle(const ForceField& ff,Atom* a1,Atom* a2,Atom* a3) //Assume the rule to be on a2 by default (No check)
{
	std::vector<double> p=ff.getAngleParamtersFailsafe(a1,a2,a3);
	if(angp) p[0]=angv;
	return p;
}
inline std::vector<double> chemtools::Rule::predictDihedral(const ForceField& ff,Atom* a1,Atom* a2,Atom* a3,Atom* a4) //Assume the rule to be on a2/a3 by default (No check)
{
	std::vector<double> p=ff.getDihedralParametersFailsafe(a1,a2,a3,a4);
	if(dihp) {p[0]=dihv; p[2]=dihm;}
	return p;
}
inline Atom* chemtools::Rule::suggest(Atom* a,const ForceField* ff) const {return ff->getAtom(randomSelect(group));}

/*std::vector<Bond> ForceField::length_parameters;
std::vector<std::pair<vector_s,vector_d>> ForceField::angle_parameters;
std::vector<std::pair<vector_s,vector_d>> ForceField::dihedral_parameters;
std::vector<std::pair<vector_s,pair_dd>> ForceField::improper_parameters;
std::vector<std::pair<pair_ss,pair_dd>> ForceField::one_four_interactions;
std::vector<std::pair<std::string,std::vector<double>>> ForceField::other_interactions;
std::vector<std::pair<std::string,double>> ForceField::one_four_charges;*/

void ForceField::loadLengths(const std::string& filename,bool erase)
{
	cout << "Loading lengths from: "<<filename<<"\n";
	//length_parameters=std::vector<std::pair<pair_ss,pair_dd>>();
	if(erase) length_parameters=std::vector<BondData>();
	//bondables=std::vector<pair_ss>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first,second;
	double b0,kb;
	while(1)
	{
		bool ffirst=false,fsecond=false;
		f>>first>>second>>b0>>kb;
		if(f.eof())
			break;
		for(Atom* a : atom_types)
		{
			if(a->toString()==first) ffirst=true;
			if(a->toString()==second) fsecond=true;
			if(ffirst && fsecond) break;
		}
		if(ffirst && fsecond) length_parameters.push_back(BondData(first,second,b0,kb));
		/*std::pair<std::string,std::string> ats(first,second);
		std::pair<double,double> vals(b0,kb);*/
		//bondables.push_back(ats);
	}
  cout << "Loaded length data\n";

  std::vector<std::string> aload;
  std::vector<std::vector<BondData>> data;
  int ind1=0,ind2=0;
  for(BondData& b : length_parameters)
  {
    ind1=-1,ind2=-1;
    for(int i=0;i<aload.size();i++)
    {
      if(aload[i]==b.atom1) ind1=i;
      if(aload[i]==b.atom2) ind2=i;
      if(ind1!=-1 && ind2!=-1) break;
    }
    //cout << b.atom1 <<" "<<b.atom2<<"\t"<<(ind1==-1)<<" "<<(ind2==-1)<<" for "<<aload.size()<<"\n";
    if(ind1==-1) {aload.push_back(b.atom1); data.push_back(std::vector<BondData>()); ind1=aload.size()-1;}
    data[ind1].push_back(b);
    if(b.atom1==b.atom2) continue;
    if(ind2==-1) {aload.push_back(b.atom2); data.push_back(std::vector<BondData>()); ind2=aload.size()-1;}
    data[ind2].push_back(b);
  }
  for(int i=0;i<aload.size();i++)
  {
    //cout <<aload[i]<<"\t";
    //for(BondData& d : data[i]) cout<<d.getSecond(aload[i])<<"\n";
    atombondmap.push_back(make_pair(aload[i],data[i]));
  }
  cout << "Loaded mapping for bondtypes\n";
	cout << "Done with lengths\n";
	f.close();
}
void ForceField::loadAngles(const std::string& filename,bool erase)
{
	cout << "Loading angles from: "<<filename<<"\n";
	if(erase) angle_parameters=std::vector<std::pair<vector_s,vector_d>>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first,second,third;
	double p1,p2,p3,p4;
	int k=0;
	bool fnd;
	while(1)
	{
		f>>first>>second>>third>>p1>>p2>>p3>>p4;
		p1*=PI/180.0; //Convert to radians
		if(f.eof())
			break;
		std::vector<std::string> ats; ats.push_back(first); ats.push_back(second); ats.push_back(third);
		std::vector<double> vals; vals.push_back(p1); vals.push_back(p2); vals.push_back(p3); vals.push_back(p4);
		angle_parameters.push_back(make_pair(ats,vals));
		first=first[0]; second=second[0]; third=third[0];
		std::vector<std::string> atg; atg.push_back(first); atg.push_back(second); atg.push_back(third);
		fnd=false;
		for(auto& a : angle_parameters)
		{
			if(matches(get<0>(a),atg)) {fnd=true; break;}
		}
		if(!fnd)
			angle_parameters.push_back(make_pair(atg,vals));
	}
	f.close();
	cout << "Done with Angles\n";
}
void ForceField::loadDihedrals(const std::string& filename,bool erase)
{
	cout << "Loading dihedrals from: "<<filename<<"\n";
	if(erase) dihedral_parameters=std::vector<std::pair<vector_s,vector_d>>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first,second,third,fourth;
	double p1,p2,p3,d;
	bool fnd;
	while(1)
	{
		f>>first>>second>>third>>fourth>>p1>>p2>>p3;
		if(f.eof())
			break;
		p1*=PI/180.0;
		std::vector<std::string> ats; ats.push_back(first); ats.push_back(second); ats.push_back(third); ats.push_back(fourth);
		std::vector<double> vals; vals.push_back(p1); vals.push_back(p2); vals.push_back(p3);
		dihedral_parameters.push_back(make_pair(ats,vals));
		fnd=false;
		first=first[0]; second=second[0]; third=third[0]; fourth=fourth[0];
		std::vector<std::string> atg; atg.push_back(first); atg.push_back(second); atg.push_back(third); atg.push_back(fourth);
		for(auto& a : dihedral_parameters)
		{
			if(matches(get<0>(a),atg)) {fnd=true; break;}
		}
		if(!fnd)
			dihedral_parameters.push_back(make_pair(atg,vals));
	}
	f.close();
	cout << "Done with dihedrals\n";
}
void ForceField::loadImpropers(const std::string& filename)
{
	cout << "Loading Improper(dihedral)s from: "<<filename<<"\n";
	improper_parameters=std::vector<std::pair<vector_s,pair_dd>>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first,second,third,fourth;
	double p1,p2;
	bool fnd;
	while(1)
	{
		f>>first>>second>>third>>fourth>>p1>>p2;
		if(f.eof())
			break;
		p1*=PI/180.0;
		std::vector<std::string> ats; ats.push_back(first); ats.push_back(second); ats.push_back(third); ats.push_back(fourth);
		improper_parameters.push_back(make_pair(ats,make_pair(p1,p2)));
		first=first[0]; second=second[0]; third=third[0]; fourth=fourth[0];
		std::vector<std::string> atg; atg.push_back(first); atg.push_back(second); atg.push_back(third); atg.push_back(fourth);
		fnd=false;
		for(auto& a : improper_parameters)
		{
			if(matches(get<0>(a),atg)) {fnd=true; break;}
		}
		if(!fnd)
			improper_parameters.push_back(make_pair(atg,make_pair(p1,p2)));
	}
	cout << "Done with impropers\n";
}
void ForceField::loadOneFourParameters(const std::string& filename,const std::string& chargefile=OtherInterFile)
{
	cout << "Loading 1-4 Potential Parameters from: "<<filename<<"\n";
	one_four_interactions=std::vector<std::pair<pair_ss,pair_dd>>();
	//one_four_charges=std::vector<std::pair<std::string,double>>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first,second;
	double p1,p2;
	bool fnd=false;
	while(1)
	{
		f>>first>>second>>p1>>p2;
		if(f.eof()) break;
		std::pair<double,double> pars(p1,p2);
		std::pair<std::string,std::string> names(first,second);
		one_four_interactions.push_back(make_pair(names,pars));
		first=first[0]; second=second[0];
		names=make_pair(first,second);
		fnd=false;
		for(auto& a : one_four_interactions)
		{
			if(matches(get<0>(a),names)) {fnd=true; break;}
		}
		if(!fnd)
			one_four_interactions.push_back(make_pair(names,pars));
	}
	/*f.open(chargefile,ios::in);
	cout << "Loading charges from: "<<chargefile<<"\n";
	double q;
	while(1)
	{
		f>>first>>q>>p1>>p2;
		if(f.eof()) break;
		one_four_charges.push_back(make_pair(first,q));
	}*/
	cout << "Done with One-Four interaction Parameters\n";
}
void ForceField::loadOtherInteractions(const std::string& filename)
{
	cout << "Loading Other Interaction (LJ) Parameters from: "<<filename<<"\n";
	other_interactions=std::vector<std::pair<std::string,std::vector<double> >>();
	ifstream f;
	f.open(filename,ios::in);
	std::string first;
	double p1,p2,p3;
	bool fnd=false;
	while(1)
	{
		f>>first>>p1>>p2>>p3;
		if(f.eof()) break;
		std::vector<double> vals; vals.push_back(p1); vals.push_back(p2); vals.push_back(p3);
		other_interactions.push_back(make_pair(first,vals));
		fnd=false;
		first=first[0];
		for(auto& a : other_interactions)
		{
			if((get<0>(a))==first) {fnd=true; break;}
		}
		if(!fnd)
			other_interactions.push_back(make_pair(first,vals));
	}
	cout << "Done with Other (LJ) Interaction parameters\n";
}
namespace forcefx
{
	static Eigen::Vector3d getAngleForce(Atom* a1,Atom* a2,Atom* a3,const std::vector<double>& ap);
	static inline Eigen::Vector3d getAngleForce(Atom* a1,Atom* a2,Atom* a3,const std::vector<double>& ap);
	static inline Eigen::Vector3d getDihedralForce(Atom* a1,Atom* a2,Atom* a3,Atom* a4,const std::vector<double>& dp);
	static inline std::vector<Eigen::Vector3d> getImproperForce(Atom* a1,Atom* a2,Atom* a3,Atom* a4,const std::pair<double,double>& ip);
}
class Topology
{
	std::vector<Atom*> atomtypes,atoms;
	std::vector<std::pair<int,int>> bonds;
	std::vector<std::pair<double,double>> bondpars,improppars;
	std::vector<std::vector<int>> angles,dihedrals,impropers;
	std::vector<std::vector<double>> anglepars,dihpars;
	const ForceField* ff=nullptr;
public:
	Topology(Molecule* generate,const ForceField& inpff);

	void describeTopology(std::ostream& os=std::cout) const
	{
		os << "Atom types:\n";
		for(Atom* a : atomtypes) os << a->toString()<<"\n";
		os << "-----------------\n";
		os << "Bonds:\n";
		for(auto& p : bonds)
		{
			const auto& lp=ff->getLengthParameters(atoms[get<0>(p)]->toString(),atoms[get<1>(p)]->toString());
			os << atoms[get<0>(p)]->toString() <<" - "<<atoms[get<1>(p)]->toString()<<".\tEquilibrium length: "<<get<0>(lp)<<"nm, force constant:"<<get<1>(lp)<<"\n";
		}
		os << "-----------------\n";
		os << "Angles:\n";
		for(auto& v : angles)
		{
			const auto& ap=ff->getAngleParamtersFailsafe(atoms[v[0]],atoms[v[1]],atoms[v[2]]);
			os << atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<"\t with equilibrium angle: "<<(ap[0]*(180/PI))<<" degrees, force constant: "<<ap[1]<<",("<<ap[2]<<","<<ap[3]<<") are the UB parameters"<<"\n";
		}
		os << "-----------------\n";
		os << "Dihedrals:\n";
		for(auto& v : dihedrals)
		{
			const auto& dp=ff->getDihedralParametersFailsafe(atoms[v[0]],atoms[v[1]],atoms[v[2]],atoms[v[3]]);
			os << atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<" - "<<atoms[v[3]]->toString()<<" with dihedral parameters: "<<dp<<"\n";
		}
		os << "-----------------\n";
		os << "Improper (dihedral)s:\n";
		for(auto& v : impropers)
		{
			const auto& dp=ff->getImproperParametersFailsafe(atoms[v[0]],atoms[v[1]],atoms[v[2]],atoms[v[3]]);
			os << atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<" - "<<atoms[v[3]]->toString()<<" with equilibrium angle: "<<get<0>(dp)<<", force constant: "<<get<1>(dp)<<" with initial value: "<<chemtools::getImproper(atoms[v[0]],atoms[v[1]],atoms[v[2]],atoms[v[3]])<<"\n";
		}
		os << "-----------------\n";
		os << "Topology ends\n";
		os << "-----------------\n";
	}

	void writeTopologyFile(std::ostream& os=std::cout,bool forcechg=false,double tchg=0,const std::string& ligname="LIG")
	{
		os << "; Ligand Topology for '"<<ligname<<"' generated by DeNovo <DVERSION> (devel)\n\n";
		os << "[defaults]\n";
		os << "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n";
		os << "1               2               yes             1.0     1.0\n\n";
		//Defaults manages the combination rules
		os << "[atomtypes]\n";
		os << ";name   bond_type     mass     charge   ptype   sigma         epsilon\n";
		for(Atom* a : atomtypes) os << a->toString()<<"\t"<<a->toString()<<"\t"<<a->seek_mass()<<"\t"<<a->seek_charge()<<"\tA\t"<<a->seek_sigma()<<"\t"<<a->seek_epsilon()<<"\n";
		os << ";End atom types\n\n";
		os << "[moleculetype]\n";
		os << ";name            nrexcl\n";
		os << ligname<<"\t 3\n\n";
		if(forcechg)
		{
			double chg=0; for(int i=0;i<atoms.size();i++) chg+=atoms[i]->seek_charge();
			#ifndef QUIET
				cout << "Net charge on system: "<<chg<<"\n";
				cout << "Setting net charge on system to "<<tchg<<" by uniform charge modification\n";
			#endif
			chg=(tchg-chg)/atoms.size();
			for(int i=0;i<atoms.size();i++) atoms[i]->setCharge(atoms[i]->seek_charge()+chg);
		}
		os << "[atoms]\n";
		os << ";   nr  type  resi  res  atom  cgnr     charge      mass       \n";
		//char res[8]; sprintf(res,"%8s",ligname);
		//std::string resi="    1 ";
		char nr[10],type[10],atom[9],charge[16],mass[14];
		char line[100];
		for(int i=0;i<atoms.size();i++)
		{
			sprintf(nr,"%d",i+1);
			//cout<<to_string(i+1)<<"\t" << nr << "\n";
			//sprintf(type,"%6s",atoms[i]->toString());
			//sprintf(atom,"%5s",atoms[i]->toString());
			sprintf(charge,"%1.6f",atoms[i]->seek_charge());
			sprintf(mass,"%1.6f",atoms[i]->seek_mass());
			sprintf(line,"%6s%6s    1 %6s%6s%6s%14s%12s\n",nr,atoms[i]->toString().c_str(),ligname.c_str(),atoms[i]->toString().c_str(),nr,charge,mass);
			os << line;
		}
		os << ";End atoms\n\n";
		os << "[bonds]\n";
		os << ";   ai     aj funct         r            k\n";
		for(int i=0;i<bonds.size();i++)
		{
			const std::pair<int,int>& p=bonds[i];
			sprintf(nr,"%d",get<0>(p)+1);
			sprintf(type,"%d",get<1>(p)+1);
			const auto& lp=bondpars[i]; //ff->getLengthParameters(atoms[get<0>(p)]->toString(),atoms[get<1>(p)]->toString());
			sprintf(charge,"%.4e",get<0>(lp));
			sprintf(mass,"%.4e",get<1>(lp));
			sprintf(line,"%6s%7s   1    %13s%13s;",nr,type,charge,mass);
			os << line <<" "<<atoms[get<0>(p)]->toString()<<" - "<<atoms[get<1>(p)]->toString()<<"\n";
		}
		os << ";End bonds\n\n";
		os << "[pairs]\n";
		os << ";   ai     aj    funct\n";
		for(int i=0;i<dihedrals.size();i++)
		{
			const std::vector<int>& v=dihedrals[i];
			sprintf(nr,"%d",v[0]+1);
			sprintf(type,"%d",v[3]+1);
			sprintf(line,"%6s%7s    1",nr,type);
			os << line << "; "<<atoms[v[0]]->toString()<<" - "<<atoms[v[3]]->toString()<<"\n";
		}
		os << "; End pairs\n\n";
		os << "[angles]\n";
		os << ";   ai     aj     ak    funct   theta         cth         rub           kub\n";
		char nr2[8],charge2[16],mass2[16];
		for(int i=0;i<angles.size();i++)
		{
			const std::vector<int>& v=angles[i];
			sprintf(nr,"%d",v[0]+1);
			sprintf(type,"%d",v[1]+1);
			sprintf(nr2,"%d",v[2]+1);
			const auto& ap=anglepars[i];// ff->getAngleParamtersFailsafe(atoms[v[0]],atoms[v[1]],atoms[v[2]]);
			sprintf(charge,"%.4e",ap[0]*(180.0/PI));
			sprintf(mass,"%.4e",ap[1]);
			sprintf(charge2,"%.4e",ap[2]);
			sprintf(mass2,"%.4e",ap[3]);
			sprintf(line,"%6s%7s%6s      5 %13s%13s%13s%13s",nr,type,nr2,charge,mass,charge2,mass2);
			os << line << "; "<<atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<"\n";
		}
		os << "; End angles\n\n";
		os << "[dihedrals] ;propers\n";
		os << ";   ai     aj     ak    al    func      theta         cth         mul\n";
		char nr3[8];
		for(int i=0;i<dihedrals.size();i++)
		{
			const std::vector<int>& v=dihedrals[i];
			sprintf(nr,"%d",v[0]+1);
			sprintf(type,"%d",v[1]+1);
			sprintf(nr2,"%d",v[2]+1);
			sprintf(nr3,"%d",v[3]+1);
			const auto& dp=dihpars[i]; //ff->getDihedralParametersFailsafe(atoms[v[0]],atoms[v[1]],atoms[v[2]],atoms[v[3]]);
			sprintf(charge,"%.4e",dp[0]*(180/PI));
			sprintf(mass,"%.4e",dp[1]);
			sprintf(charge2,"%d",(int)dp[2]);
			sprintf(line,"%6s%7s%7s%6s     9  %13s%13s%8s",nr,type,nr2,nr3,charge,mass,charge2);
			os << line << "; "<<atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<" - "<<atoms[v[3]]->toString()<<"\n";
		}
		os << "; End dihedrals\n\n";
		/*os << "[dihedrals]; impropers\n";
		os << "; i     j     k     l    func        q0    cq\n";
		for(int i=0;i<dihedrals.size();i++)
		{
			const auto& p=improppars[i];
			if(!get<1>(p)) continue;
			const std::vector<int>& v=dihedrals[i];
			sprintf(nr,"%d",v[0]+1);
			sprintf(type,"%d",v[1]+1);
			sprintf(nr2,"%d",v[2]+1);
			sprintf(nr3,"%d",v[3]+1);
			sprintf(charge,"%1.4f",get<0>(p));
			sprintf(mass,"%.4e",get<1>(p));
			sprintf(line,"%6s%6s%6s%6s   2   %10s%16s",nr,type,nr2,nr3,charge,mass);
			os << line <<"; "<<atoms[v[0]]->toString()<<" - "<<atoms[v[1]]->toString()<<" - "<<atoms[v[2]]->toString()<<" - "<<atoms[v[3]]->toString()<<"\n";
		}*/

		os << "[ system ]\n";
		os << " "<<ligname<<"\n;End System\n\n";
		os << "[molecules]\n";
		os << "; Compound\tnmols\n";
		os <<" "<<ligname<<"  \t  1\n";
		os <<"; End molecules\n";
		os <<"; ----- Topology ends here ------\n";
	}
	std::vector<Eigen::Vector3d> calculateInternalForces(Molecule* m=nullptr,bool takeimprop=true) const;

	double calculateTotalAngleEnergy(Molecule* src=nullptr) const;
	double calculateTotalDihedralEnergy(Molecule* src=nullptr) const;

};
