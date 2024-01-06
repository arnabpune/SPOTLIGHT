#ifndef HPP_MOLECULE
#define HPP_MOLECULE 1
#include <type_traits>
#include<cmath> //Can be removed (Included in maths.h)
#include<limits>
#include <map>
#include <vector> //Can be removed (Included in maths.h)
#include <csignal>
#include <atomic>
#include <mutex>
#include "maths.h"
#include "Eigen337/Eigen/Dense" //Can be removed
#include "Eigen337/Eigen/Core" //Can be removed (included in Atom.hpp)
#include "support/tensor.hpp"
#include "forcefield/forcefield.hpp"
#include "support/physics.hpp"
#ifndef SERIAL
#include <omp.h>
#endif
#ifdef REALTIME
#include <unistd.h>
#endif
/*
	*** Remaining Flags ***
	BIG_MOLECULES - Allow loading molecules with over 99999 atoms
	SERIAL - Run the program on a single thread (no MPI)
	REALTIME - Slow down the simulation and incrase frame dumping rate to visualize growth
	MACROMODE - Mode when the classes are mostly used from DNV-Macro functions
	FRAGGROW - Allow growth using fragments also
	PLACEHOLDERMOLECULES - Molecules are just placeholders for collections of atoms - Useful to perform analytics
	NOINTERFACE - No Y/N/O interface when writing out files (deprecated)
	NOZMQ - Switch off ZMQ library linkage
	NOJSON - Switch off LibJSoN Support
	MANUAL_SCORING - Use manual scoring (default is molecule's interaction energy) to pick molecules. If set, it you need to override the required_energy() manually in generation cpp file
*/
/*
	*** Removed Flags (used constants) ***
	USE_PHYSICS - Performs a weak gradient descent minimization for each molecule (Default: Yes)
	QUIET - Suppress a lot of stdout output (Default: No)
	VARIEDCOUNT - Use a varied counting mechanism (such as when counting residues instead of atoms for amino-acids). (Default: None - Count heavy atoms)
*/

//VERSION: Omega 3.0

using namespace std;

//Mutex
std::mutex myMutex1,myMutex2;


/**This exception is thrown when trying to access atom (usually by a copy of its pointer) when the Atom object is not a part of this Molecule object*/
class AtomNotInMoleculeException :  public exception {};
/**This exception is thrown when any of the start-of-generation methods (like samplegrow()) are called on a Molecule object with no seed atoms.*/
class GrowingWithoutSeedException : public exception {};
/**This exception is thrown when a seed was preloaded when no seed was expected*/
class PreloadedSeedsException : public exception {};
/** A generic I/O exception. It might be due to file errors also*/
class GeneralIOException : public exception {};
/**This exception is thrown when the system is unable to open a file (either due to lack of permissions or because the file  doesn't exist)*/
class UnableToOpenFileException: public GeneralIOException {};
/**This exception is thrown when DeNovo is requested to load files from an unknown/unsupported format. See the PDF for supported formats, and check the spellings at the function call*/
class UnknownFileFormatException : public exception {};
/**\deprecated This exception is thrown when parameters for calculating non-bonded interactions are not loaded*/
class NonbondedInteractionsNotLoadedException: public exception {};
/**\deprecated This exception is thrown when parameters when bonding data is not loaded*/
class BondParametersNotLoadedException: public exception {};
/**\deprecated This exception is thrown when requested bond between two atom types does not exist*/
class NoSuchBondException: public exception {};
/**This exception is thrown when a molecule attempts to enter trial mode again while already in it. Trail mode begins with the call of startTrials(). It ends with either rollback() or commit()*/
class AlreadyInTrialModeException: public exception {};
/**This is an internal exception which indicates to the program that the given atom of the molecule it at a "geometric end", i.e. there is no space to add any more heavy atoms to it (following proper angle and dihedral constraints). <br/>It is handled internally by the program*/
class NoPositionAvailableException: public exception {};
/**This is an internal exception which indicates to the program that the given two atoms are not connected by any path (when excluding a particular set of atoms). It is used to travese newly formed rings (to find the atoms that form the ring) <br/>It is handled internally by the program*/
class NoPathPossibleException : public exception {};
/**A special case of the NoPathPossibleException where we are looking for an augmenting path (needed for computing bond orders)*/
class NoAugmentingPathPossibleException : public NoPathPossibleException {};
/** This exception indicates that bond orders not calculated yet in the molecule.<br/>When a method throws this exception, use calculateBondOrders() to calculate before running that method.*/
class BondOrdersUncalculatedException : public exception {};
/**This exception indicates that the Molecule object has 0 heavy atoms. It is thrown when a minimum of one heavy atom is required to be present before a certain process begins*/
class EmptyMoleculeObjectException : public exception {};
/**This exception is thrown when the loaded ForceField does not have an "aromatic" Category when it is required*/
class NoAromaticCategoryException : public exception {};
/**This exception is thrown when Varied count mode is chosen, but VCLIST (list of atom types to count based on) is not set.*/
class VariedCountListNotSetException : public exception {};
/**When running by using charges calculated externally, SMILES strings are generated. If the result of the SMILES (retrieved from the network) do not match the structures, this exception is thrown.<br/> Unless thrown manually as well, it is internally handled.*/
class UnmatchingSMILESAndStructurePairException : public exception {};
//Some constants
static const int DEFAULT_MARX=5,DEFAULT_MARY=5,DEFAULT_MARZ=5,DEFAULT_RES=10,LEVELLIMIT=20,RBDESCSTEPS=7500,FULLDESCSTEPS=2000,WAIT_TIME=20;
static const double SINGLEATOMINTENERLIM=250,SINGLEATOMENERLIM=2500,RCUT=0.13,PROTCUT=0.23,ENLJRCUT=1.0,ENELRCUT=5.0,CYCCUT=0.1805,SERCUT=1,CYCTHETACUT=(12*PI)/180.0,DIHEDENERCUT=0.025,CYCLENVAR=0.4; // 0.01 nm (10pm), 1nm cutoff of SE interaction
//Typedefs
typedef tensor3D_int Grid;
//Shared variables
static double d14=0,dot=0,Ekq=0,Elj=0,RESTRAINDISTANCE=0.8,CYCPROB=0.0;
static const int ATTEMPTTRIALS=10;
static const Category* VCLIST=nullptr;
//static double common_avgenergy;
static Molecule* NULLMOL=nullptr;
//static int ZEROINT=0;

//Atom* cycat=nullptr;
class Protein;
class System;

Protein* envprot=nullptr;

/** \class System Molecule.hpp "graph/Molecule.hpp"
The System class is designed envisaging a true physical "system" with multiple interacting components.<br/>
It is this object that also stores ambient condition data (such as temperature), and the multiple components in this System are expected to "interact" with each other<br/>
From a more algorithmic perspective, all Molecule objects in a given System have non-bonding interaction energies calculated between them (during generation). The exception to this is that two macromolecules (Protein class objects), if present in the same system do not have interaction energies calculated automatically (you can specifically perform this time-taking calculation)<br/>
This algorithm assumes all atoms of a macromolecule to be stationary during generation.<br/>
See algorithm PDF for more details on generation.
*/
class System
{
  std::vector<Molecule*> mols;
  Grid* myGrid=nullptr;
	tensor3D_double* LJGrid=nullptr;
	const quickgeom::Container* restrbox=nullptr;
public:
  Eigen::Vector3d gridOffset; /**<\deprecated Grid variables are deprecated*/
  int grid_width; /**<\deprecated Grid variables are deprecated*/
	int grid_height; /**<\deprecated Grid variables are deprecated*/
	int grid_depth; /**<\deprecated Grid variables are deprecated*/
	int grid_resolution; /**<\deprecated Grid variables are deprecated*/
	const ForceField* ff; /**<The ForceField (object pointer) linked to this System. DeNovo versions upto θ2 used a direct copy of the object.<br/> Owing to the enormous number of copies created per system while loading into a queue, this was changed to a pointer*/
	double temp=298; /**<System temperature*/
public:
	/**@brief The default constructor. Construction of a System requires a ForceField object*/
	System(const ForceField& cff) :ff(&cff) {myGrid=nullptr; mols=std::vector<Molecule*>();}
	/**@brief The default copy constructor
		 @details The Molecules of this system are not explicitly cloned. Only pointers are copied.<br/>See clone(std::vector<Molecule*>) const for other options in cloning
	*/
	System(const System& sys) : ff(sys.ff)
	{
		mols=sys.mols;
		gridOffset=sys.gridOffset;
		grid_width=sys.grid_width;
		grid_height=sys.grid_height;
		grid_depth=sys.grid_depth;
		grid_resolution=sys.grid_resolution;
		myGrid=sys.myGrid;
		restrbox=sys.restrbox;
		temp=sys.temp;
	}

	/**@brief Add a molecule to the given System.
		 @details This method adds a molecule to this system. A single molecule can be part of multiple systems, but as only pointers are held in a System, changes in one reflect on the other. <br/> See the clone() methods to avoid this.
	*/
  void addMolecule(Molecule* m) {mols.push_back(m);}
  //double calculateTotalNonBondingEnergy(Atom* dummy,Molecule* ptr) const;
	//double calculateTotalNonBondingEnergy(Molecule* ptr) const;
  //std::pair<double,double> calculateDeltaEnergy(Atom* dummy,Atom* src,const Bond& b,Molecule* srcM) const;

	/**\deprecated Grid functions are deprecated*/
  inline Grid& grid() {return *myGrid;}
	/**\deprecated Grid functions are deprecated*/
  inline const Grid& grid() const {return *myGrid;}
	/**\deprecated Grid functions are deprecated*/
  void generateGrid(int,int,int,double);
	/** @brief Allows changing the ForceField linked with this system after construction.
			@details This method allows assigning a new ForceField to the System at run-time.<br/> However, it is recommended to keep using the ForceField with which the System was constructed throughout its lifetime.
	*/
	inline void setForceField(const ForceField* ffd) {ff=ffd;}
	/**@brief Set the temperature of the system. Default temperature is 298K*/
	inline void setTemperature(const double& t) {temp=t;}
	/**@brief Set the container for restraining growth to within it
		 @details The Container (supplied as a pointer) will be the region (if set) to which growth is restrained. You can further restrain it using hotspot atoms as well.
	*/
	inline void setContainer(const quickgeom::Container* cnt) {restrbox=cnt;}
	inline const quickgeom::Container* getContainer() const {return restrbox;}
	/**\deprecated Grid functions are deprecated*/
	inline bool occupied(int x,int y,int z) const {return (bool)((myGrid)?(*myGrid)(x,y,z):false);}
	/**@brief Get all the molecules contained in this System. (alternate method for a "const" System object)
		 @details This method returns all the Molecule objects (as pointers) contained in this System.<br/> Any changes made to molecules returned by this method will reflect on the molecules in the system as well.<br/> Use Molecule::Molecule(const Molecule*) to manually copy a molecule to avoid this.
		 @return std::vector of Molecule object pointers (std::vector<Molecule*>)
	*/
	inline const std::vector<Molecule*>& getMolecules() const {return mols;}
	/**@brief Get all the molecules contained in this System.
		 @details This method returns all the Molecule objects (as pointers) contained in this System.<br/> Any changes made to molecules returned by this method will reflect on the molecules in the system as well.<br/> Use Molecule::Molecule(const Molecule*) to manually copy a molecule to avoid this.
		 @return std::vector of Molecule object pointers (std::vector<Molecule*>)
	*/
	inline std::vector<Molecule*>& getMolecules() {return mols;}
	/**@brief Make a copy of this system
		 @details This method makes a copy of the given System with all the parameters. None of the Molecules in the system are copied. The pointers point to the same object.<br/>This is the same as using the assignment operator<br/> See also: clone(std::vector<Molecule*> omit) const and clone(Molecule* omit) const
	*/
	inline System clone() const;
	/**@brief Make a copy of this system. Explicitly copy some molecules as well
		 @details This method makes a copy of the given System with all the parameters. Specified molecules (if they are in the system) are also separately cloned and added to the new System. Other Molecule pointers point to their old objects.
	*/
	System clone(std::vector<Molecule*> omit) const;
	/**@brief Make a copy of this system. Explicitly copy a particular Molecule
		 @details This method makes a copy of the given System with all the parameters.<br/>It is basically a wrapper for calling clone(std::vector<Molecule*> omit) with a single Molecule pointer
	*/
	inline System clone(Molecule* omit) const {std::vector<Molecule*> tmp; tmp.push_back(omit); return clone(tmp);}
	/**\deprecated Grid functions are deprecated*/
	Eigen::Vector3i getClosestFreeSpace(int x,int y,int z,int maxtries);

	/**@brief Dump all the molecules in this System sequentially into a file
		 @param[in] output_file: Name of the file to which output must be written (type: std::string)
	*/
	void dumpSystem(const std::string& output_file) const;
	/**@brief Dump all the molecules in this System sequentially into a file. Dump TXT files (checkpoints) into a separate files.
		 @details This method creates two files. One is the normal PDB dump with multiple molecules. The other is a TXT file (the original format for DeNovo). <br/> Newer versions of DeNovo can rad directly from PDB files, and hence this method is rarely used. However, it is kept for backward support.<br/>TXT files are the standard (and only) input for all DeNovo versions before SIGMA (β,Ɣ(++), δ(++(H)), η(++))
		 @param[in] output_file: Name of the file to which normal (PDB) output must be written (type: std::string)
		 @param[in] txt_file: Name of the file to which TXT output must be written (type: std::string)
	*/
	void dumpSystem(const std::string& output_file,const std::string& txt_file) const;
	/** Get the first non-macromolecule object from the list of molecules in this System*/
	Molecule* getDrugMolecule() const;
	/** Get the first non-macromolecule object from the list of molecules in this System and remove it from the System*/
	Molecule* popDrugMolecule();
	/** Get the first macromolecule object (Protein) from the list of molecules in this System*/
	Protein* getProteinMolecule() const;
	/**@brief Empty this system (Remove all molecules)*/
	void emptySystem() {mols=std::vector<Molecule*>();}

	/** @brief Checks if a molecule is clashing with any of the molecules in the System
			@details This method uses vanderwaal radii of all the atoms to check if the Molecule being considered clashes with any other in the System. <br/> The input molecule does not necessarily have to belong to the System<br/> See also: strongWVCheck(Atom*) const
			@param[in] exm: The molecule being examined. (Molecule*)
	*/
	bool strongWVCheck(Molecule* exm) const;
	/** @brief Checks if an atom is clashing with any of the molecules in the System
			@details This method uses vanderwaal radii of all the atoms to check if the Atom being considered clashes with any other in the System. <br/> The input molecule <b>must not</b> belong to any molecule present in the System<br/> See also: strongWVCheck(Molecule*) const
			@param[in] exm: The atom being examined. (Atom*)
	*/
	bool strongWVCheck(Atom* exm) const;

	System& operator=(const System& sys)
	{
		assert(sys.ff==ff);
		mols=sys.mols;
		gridOffset=sys.gridOffset;
		grid_width=sys.grid_width;
		grid_height=sys.grid_height;
		grid_depth=sys.grid_depth;
		grid_resolution=sys.grid_resolution;
		myGrid=sys.myGrid;
		restrbox=sys.restrbox;
		temp=sys.temp;
		return *this;
	}
};

//End System

//Energy functions
/**@brief Uses σ and ε values to calculate L-J potential at a given separation 'r'*/
inline double calculateLJEnergy(const double& s1,const double& s2,const double& e1,const double& e2,const double& r,const double& dev=0)
{
  if(r>ENLJRCUT) return 0;
  //return 4*sqrt(e1*e2)*(pow(((s1+s2)/2.0)/r,12)-pow(((s1+s2)/2.0)/r,6));
	return 4*sqrt(e1*e2)*(pow(((s1+s2)/2.0)/r,12)*(1+(156*dev*dev)/(r*r))-pow(((s1+s2)/2.0)/r,6)*(1+(42*dev*dev)/(r*r)));
}
/**@brief Calculates electrostatic potential using two (partial) charges and a separation distance*/
inline double calculateElectrostaticPotential(const double& q1,const double& q2,const double& r,const double& dev=0)
{
  if(r>ENELRCUT) return 0;
	#ifndef STEP_DIELECTRIC
  	return ((DIEL*q1*q2)/r)*(1+(dev*dev)/(2*r*r));
	#else
		return ((((r<DIELST)?DIELB:DIELT)*q1*q2)/r)*(1+(dev*dev)/(2*r*r));
	#endif
}
//Force functions: Note that this force is in unconverted (raw) units. The unit obtained from here is 10^12 N/mol. The reason is that dividing with mass (in a.m.u) directly provides acceleration in nm/s^2
inline double calculateLJForce(const double& s1,const double& s2,const double& e1,const double& e2,const double& r)
{
	if(r>ENLJRCUT) return 0;
	return -(4*sqrt(e1*e2)*(pow(((s1+s2)/2.0)/r,12)*(12/r)-pow(((s1+s2)/2.0)/r,6)*(6/r)));
}
inline double calculateElectrostaticForce(const double& q1,const double& q2,const double& r)
{
	if(r>ENELRCUT) return 0;
  return -((DIEL*q1*q2)/(r*r));
}

//Extra ease of access methods
	//std::pair
/*template<class T,class U> inline U getSecond(std::pair<T,U> p,const U& u)
{
	if((get<0>(p))==u) return (get<1>(p));
	else if((get<1>(p))==u) return (get<0>(p));
	else {cout <<"Warning: None in pair (getSecond: ~Line 50 in Molecule.cpp)";}
}*/
//Template comparisons
template<class T,class U> std::pair<U,T> swap(const std::pair<T,U>& p) {return make_pair((get<1>(p)),(get<0>(p)));}

//User friendly prompts
/**@brief User-friendly prompts regarding overwriting existing files
	 @details This function takes a filename and tries to create a new file with that name. If the file exists, the user is prompted to choose one option<br/>
	 <ul>
	 	<li>'A' (Append): Append to the existing file</li>
		<li>'O' (Overwrite): Replace the file with a new one </li>
		<li>'R' (Rename): Use another name for the file</li>
		<li>'S' (Skip file generation): Skip the file opening process completely. The stream will not be functional. Note, this may lead to logging problems (needs to be fixed)
		<li>'H' (Halt): Interrupt the program for user intervention. The program will exit unconditonally, and will need to be re-run</li>
	 </ul>
*/
ofstream& fileExists(const std::string& logfname,char mode='a');

/**
	The chemtools namespace is a holder for quick functions to be used in the program (and outside if required).
	It basically contains wrapper functions to directly use Atom objects and accurately pass on their parameters to the correct candidate functions.
*/
namespace chemtools
{
	/**@brief Predict allowed positions based on angle and dihedral parameters
	 	 @details When a new atom if to be connected to an existing atom, and the trialling of the bond is being done, this function trials multiple positions, and provides a list of potential candidate posiitions, and the interaction energy at each spot.
		 @param[in] src: The source atom which is part of the newly forming bond.
		 @param[in] mol: The source molecule (which must also contain src).
		 @param[in] natom: The new atom-type being trialled.
		 @param[in] ff: The System under consideration. It is used only to get the working ForceField, and hence the name "ff".
		 @param[in] num: Number of trial positions to generate (may increase multiplicatively based on structure generated until now, and position of src in the molecule - See algorthim PDF for details on position trialling).
		 @param[in] tries: If 'num' positions could not be generated in 'tries' trials, stop.
		 @param[in] restr: Restrain the positions? (See algorithm PDF, version σ++ for more details on restrained growth).
		 @param[in] aromcat: A pointer to the aromatic category. Starting from DeNovo version θ2, the dihedral is forced to 0 if the source atom (src) is aromatic, but hasn't yet cyclized.
	*/
	static std::vector<std::pair<Eigen::Vector3d,std::pair<double,std::pair<double,double>>>> getAllAllowedPositions(Atom* src,Molecule* mol,Atom* natom,const System& ff,int num=10,int tries=10,bool restr=false,const Category* aromcat=nullptr,int devid=0);
	/**@brief Get a (boltzmann-factor weighted) random positions from all the trialled positions. See (chemtools::getAllAllowedPositions)*/
	static std::pair<Eigen::Vector3d,std::pair<double,double>> getAllowedPositions(Atom* src,Molecule* mol,Atom* natom,const System& ff,int num=10,int tries=10,bool restr=false,const Category* aromcat=nullptr,int devid=0);
	/**@brief Get the L-J (Lennard Jones) interaction energy between two atoms. It uses the internal σ, ε from the atom-type data*/
	static inline double getLJInteraction(Atom* a1,Atom* a2) {return calculateLJEnergy(a1->seek_sigma(),a2->seek_sigma(),a1->seek_epsilon(),a2->seek_epsilon(),a2->distanceFrom(a1),(a1->posdev+a2->posdev)/2);}
	/**@brief Get the L-J (Lennard Jones) force between two atoms. It uses the internal σ, ε from the atom-type data*/
	static inline double getLJForce(Atom* a1,Atom* a2) {return calculateLJForce(a1->seek_sigma(),a2->seek_sigma(),a1->seek_epsilon(),a2->seek_epsilon(),a2->distanceFrom(a1));}
	/**@brief Get the electrostatic interaction energy between two atoms. It uses the charges from the atom-type data*/
	static inline double getElectrostaticEnergy(Atom* a1,Atom* a2) {return calculateElectrostaticPotential(a1->seek_charge(),a2->seek_charge(),a2->distanceFrom(a1),(a1->posdev+a2->posdev)/2);}
	/**@brief Get the electrostatic force between two atoms. It uses the charges from the atom-type data*/
	static inline double getElectrostaticForce(Atom* a1,Atom* a2) {return calculateElectrostaticForce(a1->seek_charge(),a2->seek_charge(),a2->distanceFrom(a1));}
	/**@brief Get the electrostatic interaction energy between two atoms distance 'r' apart.
		 @details This function does not use the positions of the atoms provided to get the distance. Instead these atoms are used only as holders for atom-type data and it uses the charges from the atom-type data as well.
	*/
	static inline double getElectrostaticEnergy(Atom* a1,Atom* a2,double r) {return calculateElectrostaticPotential(a1->seek_charge(),a2->seek_charge(),r,(a1->posdev+a2->posdev)/2);}
	/**@brief Get the electrostatic force between two atoms distance 'r' apart.
		 @details This function does not use the positions of the atoms provided to get the distance. Instead these atoms are used only as holders for atom-type data and it uses the charges from the atom-type data as well.
	*/
	static inline double getElectrostaticForce(Atom* a1,Atom* a2,double r) {return calculateElectrostaticForce(a1->seek_charge(),a2->seek_charge(),r);}
	/**@brief Get the L-J (Lennard Jones) interaction energy between two atoms distance 'r' apart.
		 @details This function does not use the positions of the atoms provided to get the distance. Instead these atoms are used only as holders for atom-type data and it uses the internal σ, ε from the atom-type data as well
	*/
	static inline double getLJInteraction(Atom* a1,Atom* a2,double r) {return calculateLJEnergy(a1->seek_sigma(),a2->seek_sigma(),a1->seek_epsilon(),a2->seek_epsilon(),r,(a1->posdev+a2->posdev)/2);}
	/**@brief Get the L-J (Lennard Jones) force between two atoms distance 'r' apart.
		 @details This function does not use the positions of the atoms provided to get the distance. Instead these atoms are used only as holders for atom-type data and it uses the internal σ, ε from the atom-type data as well
	*/
	static inline double getLJForce(Atom* a1,Atom* a2,double r) {return calculateLJForce(a1->seek_sigma(),a2->seek_sigma(),a1->seek_epsilon(),a2->seek_epsilon(),r);}
	/**@brief A wrapper that uses the stored parameters of Atom objects to calculate non-bonding interaction between two specific atoms (using their stored co-ordinates for distance calculation). See Also: getNonbondingPotential(Atom*,Atom*,double)*/
	static inline double getNonbondingPotential(Atom* a1,Atom* a2) {return getLJInteraction(a1,a2)+getElectrostaticEnergy(a1,a2);}
	/**@brief A wrapper that uses the stored parameters of Atom objects to calculate non-bonding interaction between two specific atom-types at a given distance*/
	static inline double getNonbondingPotential(Atom* a1,Atom* a2,double r)
  {
    //if(r>ENRCUT) return 0;
    return getLJInteraction(a1,a2,r)+getElectrostaticEnergy(a1,a2,r);
  }

	/**@brief */
	static std::vector<std::pair<std::vector<Eigen::Vector3d>,double>> getAllAllowedFragmentPositions(Atom* src,Molecule* mol,MolecularFragment* frag,std::pair<Atom*,std::vector<Atom*>> handle,const System& ff,int num=10,int tries=10,bool restr=false,const Category* aromcat=nullptr,int devid=0);
	static std::pair<Pose,std::pair<double,double>> getAllowedFragmentPositions(Atom* src,Molecule* mol,std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> frag,const System& s,int num=20,int tries=10,bool restr=false,const Category* aromcat=nullptr,int devid=0);
}
namespace molfx
{
	#ifndef NOJSON
	void assignChargesFromJSON(Json::Value& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen=nullptr);
	#endif
	/**@brief A function that allows assignment of charges from a list obtained from an external source*/
	void assignChargesFromData(const std::vector<std::pair<std::string,double>>& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen=nullptr);
	/**@brief Same as makeCheckpoint(const std::string&,const std::vector<System>&, bool) except that it writes to any out-stream*/
	static void makeCheckpoint(std::ostream& os,const std::vector<System>& syses,bool includeprot=false);
	/**@brief Same as makeCheckpoint(const std::string&,const std::vector<System>&, bool) except that it writes to any out-stream, and uses molecule list instead of a list of systems.*/
	static void makeCheckpoint(std::ostream& os,const std::vector<Molecule*>& mols);
	/**@brief Dump a list of molecules (from a collection of "System"s) to a DeNovo checkpoint file
		 @details A DeNovo checkpoint file is an organized TXT format (See algorithm PDF), which is concatenated if multiple molecules/frames are present
	*/
	static void makeCheckpoint(const std::string& fn,const std::vector<System>& syses,bool includeprot=false)
	{
		std::ofstream fs; fs.open(fn);
		makeCheckpoint(fs,syses,includeprot);
	}
	/**@brief Load from a DeNovo checkpoint file
		 @details A DeNovo checkpoint file is an organized TXT format (See algorithm PDF), which is concatenated if multiple molecules/frames are present
	*/
	static std::vector<System> loadFromCheckpoint(const std::string& cptfile,const ForceField& ff);
	/**@brief Load multiple molecules as a list*/
	static std::vector<Molecule*> loadManyGROMolecules(std::istream& inf,const ForceField& ff,bool eval=false,bool warnf=false);
}
static ostream& operator<<(ostream& os,const chemtools::Rule& r)
{
  os << "RULE("<<r.getCentre()<<"): ";
  os << "On group: {"; for(const std::string& s : r.getRuleGroup()) os << s <<" "; os <<"}\t";
  os << "min: "<<r.getMin();
  if(r.getMax()<5) os<<", max: "<<r.getMax();
  if(r.angp) os << " with expected angle: "<<r.getAngle();
  if(r.dihp) os << " with expected dihedral: "<<r.getDihedral();
  return os;
}

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#######################Define the MOLECULE class############################|
|########################To store series of atoms############################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/


/*template<typename OT, typename RT, typename ... A>
struct lambda_expression {
    OT _object;
    RT(OT::*_function)(A...)const;

    lambda_expression(const OT & object)
        : _object(object), _function(&decltype(_object)::operator()) {}

    RT operator() (A ... args) const {
        return (_object.*_function)(args...);
    }
};*/
//Other Methods
/*inline Eigen::Vector3i castToGrid(const Eigen::Vector3d&,const System& g);
inline Eigen::Vector3i snapToGrid(const Eigen::Vector3d& r,const System& g) {return castToGrid(r,g);}
inline Eigen::Vector3d popFromGrid(const Eigen::Vector3i& r,const System& g);*/

typedef std::vector<Atom*> AtomSet;
/**Still under development. Merge two molecules into a single one*/
static Molecule* merge(const Molecule* m1,const Molecule* m2);
const Molecule* currMol=nullptr;

/** \class Molecule Molecule.hpp "graph/Molecule.hpp"
	The class is the backbone of the DeNovo program. It contains all information about the molecule and all generation methods are included here.
	 As a lot of methods in the class deal with energy, it is useful to mention that the unit of energy used in DeNovo is (kJ/mol)<br/>
	 <b>Atoms</b><br/>
	 This class holds an array (technically a 'std::vector') of Atom* (pointers to an Atom object). This storage method automatically indexes every atom by an integer.<br/>
	 When loaded from a file (PDB, GRO, or TXT), the sequence of indexing for the atoms matches the order in which they were found in the file.
	 When manually constructed, the order is same as the order in which each atom was added. If two molecules are "merged" (see merge(const Molecule*,const Molecule*)), the lists are concatenated in the same order as the parameters<br/>
	 <b>Bonds and adjacency lists</b><br/>
	 The molecule is stored in a tree format, using an adjacency list to store bonding data. Bonds are stored as vectors of integer vectors (to represent the matrix).<br/>
	 Each element (number) of an integer vector refers to the index of an Atom in the Molecule. Such vectors are compiled in another vector and indexed by the source atom (the atom which is bonded to all atoms referred to by this list)<br/>
	 A sample adjacency list (for benzene):<br/>
	 @code{.txt}
	 Atom-name,	Index,	corresponding adjacency vector

	 CA,	0,	{1,5,6}
	 CA,	1,	{0,2,7}
	 CA,	2,	{1,3,8}
	 CA,	3,	{2,4,9}
	 CA,	4,	{3,5,10}
	 CA,	5,	{4,0,11}
	 HA,	6,	{0}
	 HA,	7,	{1}
	 HA,	8,	{2}
	 HA,	9,	{3}
	 HA,	10,	{4}
	 HA,	11,	{5}
	 @endcode
	 <br/>
	 Bond orders are not stored. They usually can be determined from atom-type valencies (and sometimes with the help of hybridization)<br/>

	 <p>This class makes use of the Eigen libraries. See Eigen: http://eigen.tuxfamily.org/ for more details. (git repo: https://gitlab.com/libeigen/eigen.git)</p>
*/
class Molecule
{
protected:
	AtomSet oldatoms,oldfreeatoms; std::vector<std::vector<int>> oldbonds;
  std::vector<int> valencies;
	std::vector<std::vector<int>> bonds;
	std::vector<std::vector<unsigned short int>> orders;
	std::vector<std::vector<Atom*>> cycles,oldcycles;
public:
	std::vector<Atom*> atoms,freeatoms;
	int HC=0; /**<Hydrogen atom count*/
	int oldHC;
		int HAC=0; /**<Heavy atom count*/
		int oldHAC=0;
	#ifdef MACROMODE
		Eigen::Vector3d storedCOG,storedCOM,storedNormal;
		bool hasCOGstored=false,hasCOMstored=false;
		bool cogtype=true,comtype=true;
		Atom* headatom=nullptr; Atom* orighead=nullptr;
		int branchID=-1;
		std::vector<Atom*> headtrace,origheadtrace,branches;
		std::string cresname="DNV";
		int cresID=1;
	#endif
	int type=0; /**<Type identifier. 0 is blank (default). 1 is protein/macromolecule*/
	int ID=-1; /**<\deprecated Unique ID of this molecule*/
	bool trialmode=false;/**<Is it in "trial mode" (See flowchart in algorithm document and startTrials())*/
	double myBE=0; /**<Stored binding energy value*/
	double myHBE=0;/**<Stored binding energy value (specific to hotspot - see Protein)*/
	double errBE=0; /**<Stored binding energy error*/
	double errHBE=0; /**<Stored binding energy error (specific to hotspot - see Protein)*/
	double oldbe;
	double oldhbe;
	double olderrbe;
	double olderrhbe;
	double pE=0; /** Extra penalization (in same units as energy)*/
	#ifdef FRAGGROW
	std::vector<MolecularFragment*> componentFrags; /**<For fragment-based generation only.*/
	#endif

protected:
	/**@brief The blank Molecule constructor. It just creates a ghost molecule object with no atom data*/
  Molecule() {}
public:
	/**@brief Load a single atom (eg. when starting generation from single atom seed)*/
	Molecule(Atom* a);
	/**@brief Load a molecule by just collecting atoms.
		 @details Useful to perform rigid-molecule operations that do not involve bonding. For problems that require bond-data, see Molecule::Molecule(std::vector<Atom*>,std::vector<std::vector<int>>) */
	Molecule(std::vector<Atom*>);
	/**@brief Load a molecule from a set of atoms and some connections*/
	Molecule(std::vector<Atom*>,std::vector<std::vector<int>> b);
	/**@brief Allocate space for 'n' atoms in this Molecule object*/
	Molecule(int n);
	/**@brief The Molecule class copy constructor
		 @details This constructor is the proper way to clone a Molecule object. <br/>Unless explicitly mentioned in the documentation of a specific method <b>DO NOT</b> dereference a Molecule* before using it. It may lead to memory corruption and undefined behaviour.<br/>
		 Sample usage:
		 @code{.cpp}
		 Molecule* srcMolecule; //Molecule to be copied
		 Molecule* clonedMolecule = new Molecule(srcMolecule)
		 @endcode
	*/
	Molecule(const Molecule* o);
	/**@brief Load a molecule from TXT format (see algorithm PDF for details on TXT format)
		 @details This is the Molecule class constructor that can load molecules from TXT format file. Recent DeNovo versions support PDB and GRO formats, however, this method is kept for backward support <br/> From DeNovo \omega onwards, it is used to load checkpoints. You can supply an open input file stream.<br/>See Also: Molecule::Molecule(const char *file,bool eval=false)
		 @param[in] f: Input stream from which to load
		 @param[in] eval: Evaluate bonds? (Default: no)
		 @param[in] cpt: Is a checkpoint file? (Default: no). If yes, it will expect n more lines with bonding data.
	*/
	Molecule(std::istream& f,bool eval=false,bool cpt=false,const ForceField* ff=nullptr);
	/**@brief Load a molecule from TXT format (see algorithm PDF for details on TXT format)
		 @param[in] file: Name of the file from which to load
		 @param[in] eval: Evaluate bonds? (Default: no)
	*/
	Molecule(const char *file,bool eval=false);
	/**@brief Load Molecule from a standard format (PDB or GRO)
		 @details This is the constructor used for constructing a Molecule from a standard format file.<br/>. It requires a force-field (See ForceField) with atom names defined (to match that from the input file)
		 @param[in] infile: The input file (PDB or GRO) format
		 @param[in] format: Specify the input format. Currently "pdb" and "gro" are supported. If nothing is specified, it is taken to be "pdb"
		 @param[in] ff: The ForceField object which has the atom-types referred to in this file
		 @param[in] eval: Evaluate bonds? If set to "true", it checks for bonds (using vanderwaal radii as criteria.<br/>For macromolecules, when loading from a PDB, as only non-bonding interactions matter, it is not necessary to calculate bonds.<br/> Only when loading small molecule, this is usually needed
		 @param[in] warnf : If set to true, warning (if any) will be printed during the loading. If not, most (weak) warnings will be omitted (Default: off)
		 @param[in] chain: You can specify the chain in this (if loading a multi-chain protein). By default it is taken to be 'A'. Note: You must load multiple chains from different files with different Chain IDs
	*/
	Molecule(const std::string& infile,const ForceField& ff,const std::string& format="pdb",bool eval=false,bool warnf=false,char chain='A');
	/**@brief An alternative to Molecule::Molecule(const std::string&,const ForceField&,const std::string&,bool,bool,char) where instead of filename, a stream is provided. Useful when multiple molecules are stored in the same file. See Also(loadFromCheckpoint())*/
	Molecule(std::istream& infile,const ForceField& ff,const std::string& format="pdb",bool eval=false,bool warnf=false,char chain='A');
	/**@brief Default destructor: <b>Deletes its atom data</b> and frees up memory*/
	virtual ~Molecule();

	/**@brief Actual function called in the TXT file loading constructor. See: Molecule::Molecule(std::istream&,bool,bool)*/
	void loadMoleculeFromTXT(std::istream& f,bool eval=false,bool cpt=false,const ForceField* ff=nullptr);
	void loadFromInputStream(std::istream& infile,const ForceField& ff,const std::string& format="pdb",bool eval=false,bool warnf=false,char chain='A');

	/**@brief Get the total number of atoms (including hydrogen atoms) in this Molecule*/
	inline int getSize() const {return atoms.size();}
	/**@brief Get the total number of atoms (excluding hydrogen atoms) in this Molecule.
		 @details Hydrogen atoms are considered "light" atoms by DeNovo. An atom type is considered to be a hydrogen atom if the first character of its atom name is 'H'<br/>See Also: getSize() const
	*/
  inline int getEffectiveSize() const
	{
		if(VARIEDCOUNT) return HAC;
		else return atoms.size()-HC;
	}
	/**@brief Get the pointer to the 'n'th atom in this Molecule (alternate method for a "const" Molecule object)*/
	inline const Atom* getAtom(int n) const {return atoms[n];}
	/**@brief Get the last atom in this molecule.
		 @details It is usually the atom last added to the Molecule.<br/> It is basically a wrapper for mol->getAtom(mol->getSize()-1);
	*/
	inline const Atom* getLastAtom() const {return atoms[atoms.size()-1];}
	inline Atom* getLastAtom() {return atoms[atoms.size()-1];}
	/**@brief Get the mass of this molecule (by adding masses of individual atoms)*/
	inline double getMass() const
	{
		double m=0;
		for(int i=0;i<atoms.size();i++) m+=atoms[i]->seek_mass();
		return m;
	}

	/**@brief Get a copy of the list of atoms belonging to this Molecule
		 @details This method will give a vector of Atom* (pointers) for all the atoms in the Molecule.<br/>Any modifications made to the atoms in this list will reflect on those in Molecule as well.<br/>To avoid this copy the Molecule. See Molecule::Molecule(const Molecule*)
		 @return List (vector) of pointers to all the atoms of this molecules with same indices as the original molecule (std::vector<Atom*>)
	*/
	inline const std::vector<Atom*>& getAtoms() const {return atoms;}
	/**@brief Get the bonding data of the molecule
		 @details This method will give a copy of the adjacency list stored in this object.
		 @return The matrix is returned as a vector of integer-entry vectors (std::vector<std::vector<int>>)
	*/
	inline const std::vector<std::vector<int>>& getBonds() const {return bonds;}
	/**@brief Get the pointer to the 'n'th atom in this Molecule*/
	inline Atom* getAtom(int x) {return atoms[x];}
	/**@brief Get the index value of an atom from it's pointer
		 @return The index of the atom. If the atom is not found, the result is '-1'
	*/
	int indexOf(Atom* p) const {for(int i=0;i<atoms.size();i++) {if(atoms[i]==p) return i;} return -1;}
	#ifdef MACROMODE
	/**@brief Get the (first) index value of an atom (by name))
		 @return The index of the atom. If the atom is not found, the result is '-1'
	*/
	int indexOf(const std::string& p) const {for(int i=0;i<atoms.size();i++) {if(atoms[i]->toString()==p) return i;} return -1;}
	#endif
	/**@brief a wrapper to get the calculated binding energy for this molecule. See Protein::calculateNonBondingEnergy(Molecule*).*/
	inline double binding_energy() const {return myBE;}
	/**@brief a wrapper to get the calculated binding energy (the contribution of the "hotspot" region) for this molecule. See Protein::calculateNonBondingEnergy(Molecule*).*/
	inline double hotspot_binding_energy() const {return myHBE;}
private:
  inline void set_binding_energy(double e) {myBE=e;}
  inline void set_hotspot_binding_energy(double e) {myHBE=e;}
protected:
	/**@brief a wrapper to add some amount to the binding energy stored in this Molecule.
		 @details This method is used internally by the program to increment binding energy for every atom added during the growth phase.<br/>Do not use this function from an external source.
	*/
	inline void addBE(double x) {myBE+=x;}
	/**@brief a wrapper to add some amount to the binding energy (contributed by the "hotspot") stored in this Molecule.
		 @details This method is used internally by the program to increment binding energy for every atom added during the growth phase.<br/>Do not use this function from an external source.
	*/
	inline void addHBE(double x) {myHBE+=x;}
public:
	/**@brief Add a new atom to the Molecule.
		 @details This method allows manually adding atoms to the Molecule class.<br/> Unless specifically asked for, no calculations or validations are performed on this Atom. Also, until manual bonding definition (see addBond()) is given, the Atom is not bonded to any other in the molecule irrespective of its proximity to other atoms
	*/
  inline void addAtom(Atom* a) {atoms.push_back(a); bonds.push_back(std::vector<int>()); a->setMolecule(this);}
  friend Molecule* merge(const Molecule* m1,const Molecule* m2);
	/**@brief Calculate bond orders for the bonds of this molecules (experimental)
		 @details This method parses the adjacency list and uses the atom type definition data (such as hybridization and maximium valency) to judge the orders of various bonds.<br/>Calling this method is a prerequisite for getSMILES()
		 @param[in] ff: The ForceField used to get hybridization data and valencies for estimating bond-orders.
		 @param[in] assume_complete: This parameter (if set to "true" assumes all valencies of the Molecule to be satisfied. It is useful to use this option when working with external (non CHARMM-27) valency satisfied molecules using an "empty forcefield".<br/>
		 @param[in] assign_valence: This parameter (by default: off), if switched on, recalculates the free valencies for the atoms of this Molecule. It is useful when using the "empty" force-field on an external molecule
		 See Also: Force-Field (not the class) in PDF
	*/
	void calculateBondOrders(ForceField& ff,bool assume_complete=false,bool assign_valence=false)
	{
		if(!ff.hasCategory("aromatic")) throw NoAromaticCategoryException();
		const Category& aromcat=ff.getCategory("aromatic");
		orders=std::vector<std::vector<unsigned short int>>();
		std::vector<int> exB;
		bool rep=false;
		for(int i=0;i<atoms.size();i++) {std::vector<unsigned short int> aord; for(int j=0;j<bonds[i].size();j++) aord.push_back(0); orders.push_back(aord); exB.push_back(0);}
		for(int i=0;i<atoms.size();i++)
		{
			std::vector<Atom*> a1bonds=this->getBondedAtoms(atoms[i]);
			for(int j=0;j<a1bonds.size();j++)
			{
				if(orders[i][j]!=0 && !rep) continue;
				Atom* ca=a1bonds[j];
				int ind=indexOf(ca);
				std::vector<Atom*> a2bonds=this->getBondedAtoms(ca);
				int ord=1;
				if(assume_complete)
				{
					ord=(min(1-exB[i]+(ff.getDefaultConnectivity(atoms[i])-a1bonds.size()),1-exB[ind]+(ff.getDefaultConnectivity(ca)-a2bonds.size())));
					exB[i]+=(ord-1);
					exB[ind]+=(ord-1);
				}
				else
				{
					int exv1,exv2;
					switch(atoms[i]->toString()[0])
					{
						case 'P':
							exv1=1;
							break;
						case 'S':
							exv1=2;
							break;
						default:
							exv1=0;
					}
					switch(ca->toString()[0])
					{
						case 'P':
							exv2=1;
							break;
						case 'S':
							exv2=2;
							break;
						default:
							exv2=0;
					}
					//int exv1=(atoms[i]->toString()[0]=='S' || atoms[i]->toString()[0]=='P')?1:0,exv2=(ca->toString()[0]=='S' || ca->toString()[0]=='P')?1:0;
					cout << atoms[i]->toString() << ","<<atoms[ind]->toString() << "\t" << (1+(3-atoms[i]->getHybridization())+exv1)<<"("<<(-exB[i])<<")" << "," << (1+(3-ca->getHybridization())+exv2)<<"("<<(-exB[ind])<<")" << "\n";
					ord=min(1+(3-atoms[i]->getHybridization())+exv1-exB[i],1+(3-ca->getHybridization())+exv2-exB[ind]);
					if(!atoms[i]->getHybridization() || !atoms[ind]->getHybridization() || (aromcat.contains(atoms[i]) && atoms[ind]->getStandardValency()>1) || (aromcat.contains(atoms[ind]) && atoms[i]->getStandardValency()>1)) ord=1;
					//if(atoms[i]->getStandardValency()<ord || ca->getStandardValency()<ord) ord=min(atoms[i]->getStandardValency(),ca->getStandardValency());
					exB[i]+=(ord-1);
					exB[ind]+=(ord-1);
				}
				orders[i][j]=ord;
				for(int k=0;k<a2bonds.size();k++) {if(a2bonds[k]==atoms[i]) orders[ind][k]=ord;}
			}
			if(rep) rep=false;
			else
			{
				if(atoms[i]->getStandardValency()==1 || atoms[i]->getHybridization()==3 || aromcat.contains(atoms[i])) continue;
				int s=0;
				/*if(atoms[i]->toString()[0]=='N') s++;
				else if(atoms[i]->toString()[0]=='N') s+=2;*/
				for(int t=0;t<orders[i].size();t++) s+=orders[i][t];
				cout << atoms[i]->toString() << ":"<< s <<" "<< ff.getDefaultConnectivity(atoms[i]) << "\n";
				if(s<ff.getDefaultConnectivity(atoms[i]))
				{
					rep=true;
					for(int j=0;j<a1bonds.size();j++)
					{
						int ind=indexOf(a1bonds[j]);
						exB[ind]=0;
						for(int t=0;t<orders[ind].size();t++)
						{
							orders[ind][t]=1;
							for(int u=0;u<orders[bonds[ind][t]].size();u++) if(bonds[bonds[ind][t]][u]==ind) {orders[bonds[ind][t]][u]=1; break;}
						}
					}
					i--;
				}
			}
		}
		if(assign_valence)
		{
			int s=0;
			for(int i=0;i<orders.size();i++)
			{
				s=0;
				for(int j=0;j<orders[i].size();j++) s+=orders[i][j];
				//cout << atoms[i]->toString() << "\t"<<(s-(int)orders[i].size())<<" "<<atoms[i]->seek_valency()<<"\n";
				atoms[i]->setValency(ff.getDefaultConnectivity(atoms[i])-s);
			}
			calculateFreeAtoms();
		}
	}
	/**@brief \deprecated An update on the original calculateBondOrders() method. It is currently deprecated. Use version 3*/
	void calculateBondOrders2(ForceField& ff,bool assume_complete=false,bool assign_valence=false)
	{
		if(!ff.hasCategory("aromatic")) throw NoAromaticCategoryException();
		const Category& aromcat=ff.getCategory("aromatic");
		orders=std::vector<std::vector<unsigned short int>>();
		for(int i=0;i<atoms.size();i++)
		{
			orders.push_back(std::vector<unsigned short int>());
			for(int j : bonds[i]) orders[i].push_back(0);
		}
		bool done=false,edit=false,trigger=false,trigger2=false;
		int fails=3;
		while(!done)
		{
			done=true;
			edit=false;
			for(int i=0;i<atoms.size();i++)
			{
				for(int j=0;j<bonds[i].size();j++)
				{
					if(orders[i][j]) continue;
					done=false;
					Atom* c1=atoms[i];
					Atom* c2=atoms[bonds[i][j]];
					if(!c1->getHybridization() || !c2->getHybridization() || (c1->getHybridization()==3 && c2->getHybridization()==3) || (c1->getHybridization()==3 && c1->toString()[0]!='S' && c1->toString()[0]!='P') || (c2->getHybridization()==3 && c2->toString()[0]!='S' && c2->toString()[0]!='P'))
					{
						orders[i][j]=1; edit=true;
						for(int k=0;k<bonds[bonds[i][j]].size();k++) {if(bonds[bonds[i][j]][k]==i) {orders[bonds[i][j]][k]=orders[i][j]; break;}}
						continue;
					}
					else
					{
						cout << "Checking: "<<atoms[i]->toString()<<" with "<<c2->toString()<<"\tOrders: "<<orders[i]<<"\t";
						bool valid=true;
						int vacan=0,vacan2=0;
						for(int k=0;k<orders[i].size();k++)
						{
							if(orders[i][k]>1) valid=false;
							else if(!orders[i][k]) vacan++;
						}
						for(int k=0;k<orders[bonds[i][j]].size();k++)
						{
							if(orders[bonds[i][j]][k]>1) valid=false;
							else if(!orders[bonds[i][j]][k]) vacan2++;
						}
						cout << ((valid)?"Valid, ":"Invalid, ") << vacan<<", "<<vacan2<<"\t"<<"p1: "<<((valid && c1->getHybridization()!=3) || c1->toString()[0]=='P' || c1->toString()[0]=='S' || c2->toString()[0]=='S')<<"\n";
						if((valid && c1->getHybridization()!=3) && ((vacan==1 || vacan2==1 || trigger) && ((!aromcat.contains(c1) || c1->toString()[0]=='C') && (!aromcat.contains(c2) || c2->toString()[0]=='C')) || trigger2))
						{
							cout << atoms[i]->toString()<<" with "<<atoms[bonds[i][j]]->toString()<<" makes it here\n";
							orders[i][j]=4-c1->getHybridization(); if(orders[i][j]<=1) orders[i][j]=2; edit=true;
							for(int k=0;k<bonds[bonds[i][j]].size();k++) {if(bonds[bonds[i][j]][k]==i) {orders[bonds[i][j]][k]=orders[i][j]; break;}}
						}
						else if(!valid && (c1->toString()[0]!='S' && c2->toString()[0]!='S'))
						{
							orders[i][j]=1; edit=true;
							for(int k=0;k<bonds[bonds[i][j]].size();k++) {if(bonds[bonds[i][j]][k]==i) {orders[bonds[i][j]][k]=orders[i][j]; break;}}
							continue;
						}
						else if((vacan==1 && c2->toString()[0]=='S') || (vacan2==1 && c1->toString()[0]=='S'))
						{
							if(c1->toString()[0]=='S' && c1->getStandardValency()==4) orders[i][j]=4-c2->getHybridization();
							else if(c2->toString()[0]=='S' && c2->getStandardValency()==4) orders[i][j]=4-c1->getHybridization();
							else orders[i][j]=1;
							for(int k=0;k<bonds[bonds[i][j]].size();k++) {if(bonds[bonds[i][j]][k]==i) {orders[bonds[i][j]][k]=orders[i][j]; break;}}
							continue;
						}
					}
				}
			}
			if(!edit)
			{
				if(!trigger) trigger=true;
				else
				{
					if(trigger2)
					{
						fails--;
						if(!fails) break;
					}
					trigger2=true;
				}
			}
			else {fails=3; trigger2=trigger=false;}
		}
		cout << "Final orders:\n";
		for(int i=0;i<atoms.size();i++)
		{
			cout <<atoms[i]->toString()<<"("<<i<<")"":\t";
			for(int j=0;j<orders[i].size();j++) cout << atoms[bonds[i][j]]->toString()<<"("<<bonds[i][j]<<")"<<"-"<<orders[i][j]<<" ";
			cout << "\n";
		}
		cout << "Ended order data\n";
	}

	std::vector<int> findAugmentedPath(const std::vector<int>& nunsat,int start,std::vector<int> path=std::vector<int>(),int lastidx=-1)
	{
		// Uses bond orders from object variable 'orders'
		path.push_back(start);
		std::cout << "PATH\t\t"<<path<<"\n";
		if(lastidx!=-1 && nunsat[start]) {if(path.size()%2==0) return path; else throw NoAugmentingPathPossibleException();}
		for(int bi=0;bi<bonds[start].size();bi++)
		{
			int ai=bonds[start][bi];
			if(!atoms[ai]->genpurpflag || contains(path,ai)) continue;
			if(nunsat[start]==0 && nunsat[ai]==0 && path.size()%2==0)
			{
				bool failr=false;
				for(int mbx=0;mbx<bonds[start].size();mbx++)
				{
					if(bonds[start][mbx]==ai)
					{
						if(orders[start][mbx]<=1)
						{
							std::cout << "PATHFAIL\t\t"<<start<<"-"<<bonds[start][mbx]<<"\n";
							failr=true;
						}
						break;
					}
				}
				if(failr) continue;
			}
			try {return findAugmentedPath(nunsat,ai,path,start);} catch(NoAugmentingPathPossibleException& ex) {}
		}
		if(path.size()%2==0 && nunsat[path[path.size()-1]]>0) return path;
		else throw NoAugmentingPathPossibleException();
	}

	/**@brief An update on the calculateBondOrders2() method. Returns false if some matches were not possible. Uses the augmenting path algorithm.*/
	bool calculateBondOrders3(ForceField& ff,bool assume_complete=true)
	{
		if(!ff.hasCategory("aromatic")) throw NoAromaticCategoryException();
		const Category& aromcat=ff.getCategory("aromatic");
		orders=std::vector<std::vector<unsigned short int>>();
		std::vector<int> nunsat,unsatids;
		for(int i=0;i<atoms.size();i++)
		{
			int nu=chemtools::getNumUnsat(atoms[i],&ff);
			nunsat.push_back(nu);
			if(nu>0) {unsatids.push_back(i); atoms[i]->genpurpflag=true;}
			else atoms[i]->genpurpflag=false;
			orders.push_back(std::vector<unsigned short int>());
			for(int bi : bonds[i]) orders[i].push_back(1);
		}
		bool allmatched=true;
		bool done=false;
		while(!done)
		{
			done=true;
			std::cout << "Currently at: "<<nunsat<<"\n";
			std::cout << "With orders : "<<orders<< "\n";
			for(int midx : unsatids)
			{
				if(nunsat[midx]==0) continue;
				std::vector<int> ap;
				std::cout << "Trying: "<<midx<<"\n";
				try {ap=findAugmentedPath(nunsat,midx);} catch(NoAugmentingPathPossibleException& ex) {continue;}
				if(ap.size()<=1) continue;
				else done=false;
				std::cout << "Got:\t"<<ap<<"\n";
				for(int i=0;i<ap.size()-1;i++)
				{
					std::cout <<"\t"<<ap[i]<<"\n";
					for(int bidx=0;bidx<bonds[ap[i]].size();bidx++)
					{
						if(bonds[ap[i]][bidx]==ap[i+1])
						{
							if(i%2==0) {orders[ap[i]][bidx]++; nunsat[ap[i]]--; nunsat[ap[i+1]]--;}
							else {orders[ap[i]][bidx]--; nunsat[ap[i]]++; nunsat[ap[i+1]]++;}
							break;
						}
					}
					for(int bidx=0;bidx<bonds[ap[i+1]].size();bidx++)
					{
						if(bonds[ap[i+1]][bidx]==ap[i])
						{
							if(i%2==0) orders[ap[i+1]][bidx]++;
							else orders[ap[i+1]][bidx]--;
						}
					}
				}
				break;
			}
		}
		for(int fv : nunsat) {if(fv!=0) return false;}
		return true;
	}

	/**@brief Mark the atoms which are cyclized.
		 @details This sets the "cyclized" flag for the atoms which have formed cycles to "True"
	*/
	void markCyclizedAtoms()
	{
		std::vector<int> traced;
		int stI=0;
		while(bonds[stI].size()<=1)
		{
			stI++;
			if(stI>=bonds.size()) break;
		}
		if(stI>=bonds.size()) return;
		std::vector<int> travpath;
		traverseMarkCyclicAtoms(stI,travpath);
	}
	/**@brief Calculate the number of rotatable bonds
		 @details Calculates N_{ROT}, or the number of rotational degrees of freedom. This is used to empirically correct for some of the entropy caused due to ligand size<br/>
		 <b>Note: </b> This assumes that the cycles have been already marked in this molecule (See: Molecule::markCyclizedAtoms())
	*/
	int calculateNRot(const ForceField* ff=nullptr,bool mod=false)
	{
		const Category* cat=(ff && ff->hasCategory("amidic"))?&(ff->getCategory("amidic")):nullptr;
		//if(cat) cout << "Found category: 'amidic'\n";
		int nrot=0;
		for(int i=0;i<atoms.size();i++)
		{
			if(atoms[i]->getStandardValency()<2 || atoms[i]->getHybridization()<2) continue; //Consider changing valency limit to 3 (Does C-O bond in R-C(R')-OH actually have a degree of freedom)
			for(int j=0;j<bonds[i].size();j++)
			{
				if(atoms[bonds[i][j]]->getStandardValency()<2 || atoms[i]->getHybridization()<2 || (mod && atoms[i]->getHybridization()<3 && atoms[bonds[i][j]]->getHybridization()<3) || (cat && ((cat->contains(atoms[i]) && atoms[bonds[i][j]]->toString()[0]=='N') || (cat->contains(atoms[bonds[i][j]]) && atoms[i]->toString()[0]=='N')) ) ) continue;
				bool all1=true,all2=true;
				for(const auto& v : cycles)
				{
					all1=true; all2=true;
					for(const Atom* ta : v)
					{
						if(ta==atoms[i]) all1=false;
						else if(ta==atoms[bonds[i][j]]) all2=false;
					}
					if(!all1 && !all2) {/*cout << "\t" << atoms[i]->toString() << " - " << atoms[bonds[i][j]]->toString()<<"\n";*/ break;}
				}
				if(!all1 && !all2) continue;
				//cout << atoms[i]->toString() << " - " << atoms[bonds[i][j]]->toString()<<"\n";
				nrot++;
			}
		}

		return nrot/2;
	}

	/**@brief Get Minimum coordinate of the protein*/
	double getMin(int idx)
	{
		double minx=std::numeric_limits<double>::max();
		for(Atom* a : this->getAtoms())
			if(a->getPosition()(idx)<minx) minx=a->getPosition()(idx);
		return minx;
	}

	/**@brief Get Maximum coordinate of the protein*/
	double getMax(int idx)
	{
		double maxx=std::numeric_limits<double>::min();
		for(Atom* a : this->getAtoms())
			if(a->getPosition()(idx)>maxx) maxx=a->getPosition()(idx);
		return maxx;
	}
private:
	void traverseMarkCyclicAtoms(int ind,std::vector<int>& visited,std::vector<int> path=std::vector<int>(),int skip=-1)
	{
		path.push_back(ind);
		visited.push_back(ind);
		std::vector<int> nTs=bonds[ind];
		if(nTs.size()<=1) return;
		for(int n : nTs)
		{
			if(n==skip) continue;
			if(contains(path,n))
			{
				int i;
				std::vector<Atom*> psh;
				for(i=path.size()-1;i>=0;i--)
				{
					psh.push_back(atoms[path[i]]);
					atoms[path[i]]->setCyclized(true);
					if(path[i]==n) break;
				}
				cycles.push_back(psh);
				std::vector<int> npath; for(int j=0;j<i;j++) npath.push_back(path[j]);
				path=npath;
				continue;
			}
			if(contains(visited,n)) continue;
			else traverseMarkCyclicAtoms(n,visited,path,ind);
		}
	}
public:
	/**@brief Print out the adjacency list to standard output*/
	void printStructureTree()
	{
		for(int i=0;i<bonds.size();i++)
		{
			cout << atoms[i]->toString()<<":\t";
			for(int j=0;j<bonds[i].size();j++)
			{
				cout <<atoms[bonds[i][j]]->toString()<<"("<<orders[i][j]<<"), ";
			}
			cout << "\n";
		}
	}

	/*void scamGasteigerConvergence(ForceField& ff,double cutoff=0.01,int convlimit=10,double total_charge=0)
	{
		std::vector<double> oldcharges;
		for(Atom* a : atoms) oldcharges.push_back(0);
		for(int i=0;i<convlimit;i++)
		{
			std::vector<double> charges=oldcharges;
			for(int i=0;i<atoms.size();i++)
			{
				Atom* a=atoms[i];
				std::vector<Atom*> bats=this->getBondedAtoms(a);
				for(int j=0;j<bats.size();j++)
				{
					int ind=indexOf(bats[j]);
					charges[i]-=ff.getScamGasteigerParameter(a)*(1+oldcharges[i]);
					charges[ind]-=ff.getScamGasteigerParameter(bats[j])*(1+oldcharges[ind]);
				}
			}
			oldcharges=charges;
		}
		for(int i=0;i<atoms.size();i++) atoms[i]->setCharge(oldcharges[i]);
	}*/
	/**@brief "Caps" all the atoms of this molecule
		 @details This method forcibly sets the free-valency count of each atom to 0, rejecting any remaining valencies. This process is referred to (here and later) as "capping"<br/>It is useful when using a seed molecule which as a lot of free valecies, but the growth is to be targetted from a particular atom<br/>See also: Atom::seal(), Molecule::free(Atom*,int)
	*/
	inline void sealAll()
	{
		for(Atom* a : freeatoms) a->seal();
		freeatoms=std::vector<Atom*>();
	}
	/**@brief Manually "frees" a capped or completely valency satisfied atom
		 @details Grants some free valencies to a given atom. No checks are performed at this stage to ensure that it does not cross the upper limit of allowed bonds for that atom.
		 @param[in] a: A pointer to the Atom whose valency is to be set. If it doesn't belong to this Molecule, AtomNotInMoleculeException is thrown
		 @param[in] val: The number of valecies to set for this atom. If it is negative, the free valecies of the target atom are set to 1 (if it was initally 0), otherwise it is unchanged. It should not be 0.
	*/
	inline void free(Atom* a,int val=-1)
	{
		assert(val); //val!=0
		if(!contains(atoms,a)) throw AtomNotInMoleculeException();
		if(!contains(freeatoms,a)) freeatoms.push_back(a);
		if(val<0) val=a->seek_valency();
		if(val==0) val=1;
		a->setValency(val);
	}
	/**@brief Manually parse the atom-list to add/remove atoms from the "free atoms" list. This is useful only when the Molecule has been manually constructed (using functions like addAtom(), makeBond(), and bondAtom())*/
	void calculateFreeAtoms()
	{
		freeatoms=std::vector<Atom*>();
		for(int i=0;i<atoms.size();i++) {if(atoms[i]->seek_valency()) freeatoms.push_back(atoms[i]);}
	}

	//Grid
	/**\deprecated Grid functions are deprecated*/
	virtual void loadGrid(System& s) {s.grid().fill_with_zero();}

	//Back to atoms
	/**@brief Get all atoms bonded to a given atom
		 @details Gives a list of pointers to all the atoms that are bonded to this atom.
		 @param [in] src: Source atom for which bonded atom list is required. If the atom is itself not in this molecule, AtomNotInMoleculeException is thrown
		 @return List of pointers to all the atoms (std::vector<Atom*>)
	*/
	const std::vector<Atom*> getBondedAtoms(Atom* src) const
	{
		std::vector<Atom*> res;
    //cout << src->toString()<<": getBondedAtoms() was called on this\n";
		int ind=indexOf(src);
    if(ind==-1)
		{
			cout << src->toString() << "("<<(long)src<<") not in: "; for(int i=0;i<atoms.size();i++) cout << atoms[i]->toString()<<"("<<(long)atoms[i]<<"),"; cout << "\n"; throw AtomNotInMoleculeException();
		}
		for(int i : bonds[ind]) res.push_back(atoms[i]);
		return res;
	}
	/*std::vector<Atom*> getBondedAtoms(Atom* src)
	{
		std::vector<Atom*> res;
		int ind=indexOf(src);
    if(ind==-1) throw AtomNotInMoleculeException();
		for(int i : bonds[ind]) res.push_back(atoms[i]);
		return res;
	}*/

	/**@brief Start trial mode.
		 @details This method marks the beginning of trial mode. It backs-up the current state of the Molecule object, and all modifications made from this point are subject to be lost (see rollback())<br/>Trial mode is actively used during generation to allow recoiling when the growth reaches a stage from where it cannot proceed. (See algorithm PDF for more details)<br/>You cannot use this method twice in a row. This method must be paired with either commit() or rollback(), and only after either of these is executed can this function be re-used. Not doing so may lead to an AlreadyInTrialModeException being thrown
	*/
	inline void startTrials(int devid=0);
	
	/**@brief Restore original molecule. All changes from the last call of startTrials() will be deleted
		 @details This method marks the end of "trial mode". All changes are reverted back and trial mode is switched off.<br/>See Also: commit()
	*/
	inline void rollback(const System& s,int devid=0);
	/**@brief \deprecated This method is not safe to use.
		 @details \deprecated Rollback partially (only remove those atom types which are unsatisfied).
	*/
  void miniRollback(const ForceField& ff)
  {
    cout << atoms.size() << " in the beginning\n";
    for(Atom* a : oldatoms)
		{
			cout << a->toString() << ","; cout << "\n";
		}
    std::vector<Atom*> newatoms;
    std::vector<std::vector<int>> newbonds;
    std::vector<std::pair<Atom*,std::vector<Atom*>>> bonddata;
    for(Atom* a : atoms) {bonddata.push_back(make_pair(a,getBondedAtoms(a)));}
    cout << "Bond data loaded\n";
    std::vector<Atom*> waste;
    bool passed=false;
    while(true)
    {
      cout << "Loop: "<<atoms.size() << " atoms\n";
      newatoms=std::vector<Atom*>();
      newbonds=std::vector<std::vector<int>>();
      for(int i=0;i<atoms.size();i++)
      {
        if(contains(oldatoms,atoms[i]) || ff.isSatisfied(atoms[i],getBondedAtoms(atoms[i])))
        {
          newatoms.push_back(atoms[i]);
          newbonds.push_back(std::vector<int>());
        }
        else
        {
          waste.push_back(atoms[i]);
          cout << atoms[i]->toString()<<" connected to:\t";
          for(int j=0;j<bonds[i].size();j++)
					{
						cout <<atoms[bonds[i][j]]->toString();cout<<"\n";
					}
        }
      }
      if(newatoms.size()==atoms.size()) break;
      cout << "\\Loop: "<<newatoms.size()<< " atoms\n";
      int K=0;
      for(int i=0;i<atoms.size();i++)
      {
        if(!contains(newatoms,atoms[i])) continue;
        for(int j : bonds[i]) {if(contains(newatoms,atoms[j])) newbonds[K].push_back(j);}
        K++;
      }
      atoms=newatoms;
      bonds=newbonds;
    }
    for(Atom* a : waste)
    {
      //if(a) cout << "Deleting: "<< a->toString()<<"\n";
      delete a;
    }
    cout << "Atoms fixed \tand\t";
    /*for(int i=0;i<atoms.size();i++) newbonds.push_back(std::vector<int>());
    for(int i=0;i<atoms.size();i++)
    {
      for(int j=i+1;j<atoms.size();j++)
      {
        if(contains(findByKey(bonddata,atoms[i]),atoms[j])) newbonds[i].push_back(j);
        if(contains(findByKey(bonddata,atoms[j]),atoms[i])) newbonds[j].push_back(i);
      }
    }
    bonds=newbonds;*/
    cout << "Bonds assigned\n";
    for(int i=0;i<atoms.size();i++)
    {
      if(contains(oldatoms,atoms[i])) continue;
      if(!bonds[i].size())
      {
        delete atoms[i];
        atoms.erase(atoms.begin()+i);
        bonds.erase(bonds.begin()+i);
        i--;
      }
    }
    for(int i=0;i<atoms.size();i++) atoms[i]->setValency(atoms[i]->getStandardValency()-bonds[i].size());
    cout << "Valencies set\n";
    freeatoms=std::vector<Atom*>();
    HC=0;
		if(VARIEDCOUNT) HAC=0;
    for(Atom* a : atoms)
    {
      if(a->canBond()) freeatoms.push_back(a);
      if(a->isHydrogen()) HC++;
			if(VARIEDCOUNT && VCLIST->contains(a)) HAC++;
    }
    cout << atoms.size() << " atoms left\n";
    //for(Atom* a : atoms) cout << a->toString() << ","; cout << "\n";
    trialmode=false;
  }
	/**@brief Save all changes and end trial mode
		 @details This method marks the end of "trial mode". All changes made are kept, the back-up copy is deleted, and trial mode is switched off.
	*/
	inline void commit();

	//Simple random method
	/**@brief Randomly select an atom with free valencies
		 @details Randomly picks an atom from the list of "free atoms". If the Molecule was manually constructed, you may need to use calculateFreeAtoms()
		 @return A pointer to a random Atom belonging to this molecule which has at-least one free valency. If there is no such atom, a nullptr (Null-pointer) is returned
	*/
	inline Atom* getRandomBondableAtom() const
	{
		if(!freeatoms.size()) return nullptr;
		return randomSelect(freeatoms);
	}

	/**@brief Join two atoms to make a bond
		 @details This is one of the functions used for manually constructing a molecule in cpp. I can also  be used to add non-existing/un-expected bonds to a loaded molecule.<br/>This function accepts two atoms (as pointers). These atoms <b>must belong to this Molecule object</b>. See addAtom(Atom*) for details on adding an atom<br/>It is necessary for both atoms to have one free valency at-least. Make sure you have chosen the correct atom types before bonding. You can add valencies using free() or Atom::setValency(int), but it is not recommended. If you do not know what atom-type to use and/or do not intend to use it with force-fields, use the "empty forcefield".<br/> See Also: Force-Field in PDF, and the ForceField class
	*/
	void makeBond(Atom* a1,Atom* a2)
	{
		assert(a1!=a2);
		int i1=indexOf(a1),i2=indexOf(a2);
		if(i1==-1 || i2==-1) throw AtomNotInMoleculeException();
		#ifndef PLACEHOLDERMOLECULES
			if(!a1->canBond() || !a2->canBond()) throw NoValencyAvailableException();
		#endif
		bonds[i1].push_back(i2);
		bonds[i2].push_back(i1);
		a1->addBond(); a2->addBond();
		if(!a1->canBond()) freeatoms.erase(std::remove(freeatoms.begin(), freeatoms.end(), a1), freeatoms.end());
		if(!a2->canBond()) freeatoms.erase(std::remove(freeatoms.begin(), freeatoms.end(), a2), freeatoms.end());
	}

	/**@brief Add an atom to the molecule by bonding it to a pre-existing atom of the molecule
		 @details This method is a short way to add an atom to a Molecule when manually building a molecule. It adds the Atom to this Molecule and then bonds it to the chosen atom (which must be in this molecule).<br/> It achieves both the steps in one function.(see addAtom(Atom*) and makeBond(Atom*,Atom*))<br/>No checks are performed to see if the newly added atom already belongs to the molecule. The use must manually check this, or this may lead to memory corruption and/or infinite loops (see Molecule::indexOf(Atom*) const).
		 @param[in] natm: The newly added atom (this must not belong to the molecule object when this method is called)
		 @param[in] src: The source atom (to which this new atom is attached. This atom must be in the Molecule. No checks are performed here either)
		 @return true/false (bool): If either atom (new or source) cannot bond due to lack of valences, no bond is formed and 'false' is returned. Otherwise the function exists with a return value 'true'
	*/
	bool bondAtom(Atom* natm,Atom* src)
	{
    //Atom* nacpy=natm;
		if(!natm->canBond() || !src->canBond()) return false;
    if(natm->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(natm)) HAC++;
		int sI=indexOf(src);
		atoms.push_back(natm);
		bonds.push_back(std::vector<int>(1,sI));
		bonds[sI].push_back(atoms.size()-1);
		src->addBond(); natm->addBond();
		if(natm->canBond()) freeatoms.push_back(natm);
		if(!src->canBond()) freeatoms.erase(std::remove(freeatoms.begin(), freeatoms.end(), src), freeatoms.end());
		return true;
	}

	/**@brief Bond a fragment to a given molecule by attaching one point to another atom. Returns a bool value to report success*/
	std::vector<Atom*> bondFragment(MolecularFragment* nfrag,int ind,Atom* src);

	/**@brief A wrapper to iterate through all the atoms in this molecule, and sum up the non-bonding interaction energy of this new (external) atom with this molecule. See Also: chemtools::getNonbondingPotential()*/
	inline virtual std::pair<double,double> calculateNonBondingEnergy(Atom* dummy,bool write=false);

	/**@brief \deprecated An old trial method for growth*/
	void samplegrowimpartial(const System& s,int n,int trials=15)
	{
		int t=0;
		while(getEffectiveSize()<n)
		{
			if(t>trials)
			{
				cout << "Too many failures: \n";
				break;
			}
			Atom* gP=getRandomBondableAtom();
			auto nbs=s.ff->getAllowedBonds(gP);
			BondData nb=randomSelect(nbs);
			Atom* nat=nb.getSecond(gP,*s.ff);
			//std::vector<Eigen::Vector3d> alposes=
			//if(!alposes.size()) {t++; continue;}
			Eigen::Vector3d pos;
      try{pos=get<0>(chemtools::getAllowedPositions(gP,this,nat,s,15,5)); t=0;}
      catch(NoPositionAvailableException& ex) {t++; continue;}
			nat->setPosition(pos);
			bondAtom(nat,gP);
			for(Atom* a : atoms)
			{
				cout << a->toString() << "("<<a->seek_valency()<<")";
			}
			cout << "\t";
			for(Atom* a : freeatoms)
			{
				cout << a->toString() << "("<<a->seek_valency()<<")";
			}
			cout <<"\n";
		}
  }
	/**@brief \deprecated Calculate the non-bonding interaction energy of the Molecule to the hotspot region of the Proteins in this system. Use calculateNonBondingEnergy()*/
  inline double getHotspotEnergy(const System& s,double v=0);
	/**@brief Attempts to add a new atom to a source atom. See flowchart of algorithm in PDF for more details*/
  bool attempt(const System& s,Atom* gP,Atom* nat,Atom*& nxt,int trials=10,bool restr=false,int nps=0,int devid=0);

	#ifdef FRAGGROW
	/**@brief Attempts to add a new atom to a source atom. See flowchart of algorithm in PDF for more details*/
  bool attemptFragment(const System& s,Atom* gP,std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> nat,Atom*& pivat,int trials=10,bool restr=false,int nps=0,int devid=0);
	#endif
	//bool join(const System& s,Atom* gP,Atom* nat,Atom*& nxt,bool acc,const std::pair<Eigen::Vector3d,std::pair<double,double>>& penv);

	/**@brief Complete all valencies by adding hydrogens
		 @details This functions adds hydrogen atoms to all the available valencies of each atom in the molecule. Geometric checks are mild and in the interest of time, exact position of hydrogen atoms are not calculated, which may lead to two overlapping hydrogen atoms from the same parent atom.<br/>The overlaps can easily be addressed by a simple steepest descent optimization in any molecule viewer (such as Chimera)<br/>The energy contribution of non-polar hydrogen atoms is very low. If you think certain hydrogen atoms have significant contribution, ensure that they are properly defined in the 'definitions' file. See algorithm PDF for more details.
	*/
  bool reduce(const System& s,int trials=15,bool forced=false,int devid=0); //Acceptable

	/**@brief \deprecated An old trial method for growth*/
  void samplegrowunrigid(const System& s,int n,int trials=15)
  {
  	int t=0;
    startTrials();
    cout << "Samplegrow begins\n";
  	while(getEffectiveSize()<n)
  	{
  		if(t>trials)
			{
				cout << "Too many failures: \n";
				break;
			}
  		Atom* gP=getRandomBondableAtom();
      if(!gP)
			{
				cout << "Out of valencies\n";
				break;
			}
      cout << gP->toString()<<"("<<gP->seek_valency()<<") chosen\n";
  		Atom* nat=s.ff->selectAtomByRule(gP,getBondedAtoms(gP),nullptr,0);
  		//if(!alposes.size()) {t++; continue;}
      Atom* nxt;
      if(!attempt(s,gP,nat,nxt)) {t++; continue;}
      bondAtom(nat,gP);
      if(nxt)
      {
        cout << nat->toString()<< "-"<<nxt->toString()<<"\n";
        makeBond(nxt,nat);
      }
      t=0;
  		for(Atom* a : atoms)
			{
				cout << a->toString() << "("<<a->seek_valency()<<")";
			}
  		cout << "\t";
  		for(Atom* a : freeatoms)
			{
				cout << a->toString() << "("<<a->seek_valency()<<")";
			}
  		cout <<"\n";
  	}
    ofstream prf; prf.open("prereduction.pdb",ios::out);
    dumpMol(prf);
    prf.close();
    reduce(s);
    ofstream rf; rf.open("postreduction.pdb",ios::out);
    dumpMol(rf);
    rf.close();
    miniRollback(*s.ff);
	}
	/**@brief Try to satisfy all atom-type requirements. See algorithm details in the PDF containing the flowchart*/
  bool completeAttempt(ostream& logfile,const System& s,Atom* src,int trials=25,int nupper=90,int level=0,int levellim=LEVELLIMIT,bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,Molecule*& common_newroot=NULLMOL,int common_seedlength=0,int devid=0);
  	/**@brief The internal growth function which guided the molecule growth
		 @details This function takes in all the necessary parameters to grow one complete molecule. It is this method which is called multiple times by other functions such as generateLigands(), quickPrune(), etc. to complete one generation.<br/>
		 It grows one molecule to the given receptor (which must be part of the System), and also stores the interaction energy of the molecule in it. <br/>
		 Also, it is this function which starts (and ends) the trial mode.
		 Refer to the algorithm flowchart for complete details of the growth algorithm.
		 @param[in] logfile: The output stream to write the log to
		 @param[in] s: The System within which this generation is placed. This provides the temperature for growth, and must also constain this molecule and the receptor.
		 @param[in] n: The The minimum size upto which to continue growth.
		 <br/>For other parameters, refer Protein::generateLigands(const ForceField&,int,int,const std::string&,std::vector<Molecule*>,double,double,int,int,int,double,bool,bool,bool,bool,Molecule*,bool,int)
	*/
  void samplegrow(ostream& logfile,const System& s,int n,int trials=15,int nupper=90,int levellim=LEVELLIMIT,bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,Molecule*& common_newroot=NULLMOL,int common_seedlength=0,int devid=0)
  {
		#ifdef REALTIME
			this->dumpMol("seed.pdb",s.ff);
		#endif
    int t=0;
    cout << "samplegrow(rigid) begins\n";
    logfile << "samplegrow(rigid) begins\n";
    bool acc=true;
    Atom* gP=nullptr;
    while(getEffectiveSize()<n)
    {
      acc=true;
      startTrials(devid);
  		if(t>trials)
			{
				cout << "Too many failures: \n"; logfile << "Too many failures: \n";
				break;
			}
			if(superseed && common_newroot && freeatoms.size() && this->getEffectiveSize()==common_newroot->getEffectiveSize()) gP=this->freeatoms[freeatoms.size()-1];
  		else gP=getRandomBondableAtom();
      if(!gP)
			{
				cout << "Out of valencies\n";
				logfile << "Ended: Out of valencies\n";
				commit();
				break;
			}
      cout << "New begin: "; describeStructure();
      if(completeAttempt(logfile,s,gP,ATTEMPTTRIALS,nupper,0,levellim,restr,superseed,reseed,deseed,common_newroot,common_seedlength,devid))
			{
				cout << "Commiting: "<<gP->toString()<<"\n";
				commit(); t=0;
			}
      else
      {
        cout << "Rollback\n";
        //this->dumpMol("pre_rolledback.pdb");
        rollback(s,devid);
				#ifdef REALTIME
        	this->dumpMol("rolledback.pdb",s.ff);
					this->dumpMol("realtime.pdb",s.ff);
					usleep(REALTIME*2);
				#endif
        t++;
        cout << "Rolled back: ";
        logfile << "Rolled back: ";
        this->describeStructure(logfile);
      }
    }
		if(this->getEffectiveSize()<n) return; //Don't waste optimization time on "small" molecules which will be rejected mostly anyway.
    //this->dumpMol("pre_reduction.pdb");
		if(USE_PHYSICS)
		{
			physics::optimizeMolecule(this,*s.ff,250);
			cout << "Short Geomtery-optimized molecule!\n";
		}
    try
		{
			this->reduce(s,trials,true,devid);
			if(FORCE_CHARGES) this->uniformChargeAdjustments(this->getNetCharge(*s.ff));
		}
    catch(AtomNotInMoleculeException& ex) {this->dumpMol("pre_reduction.pdb"); std::cout << "Crashed with AtomNotInMoleculeException\n"; throw std::exception();}
		if(USE_PHYSICS)
		{
		//physics::gradientDescentPhysics(this);
			physics::optimizeMolecule(this,*s.ff);
			cout << "Well Geomtery-optimized molecule!\n";
		}
  }

	#ifdef FRAGGROW
	MolecularFragment* completeAttemptFragment(ostream& logfile,const System& s,Atom* src,const FragmentSet& fs,int trials=25,int nupper=90,int level=0,int levellim=LEVELLIMIT,bool restr=false,int devid=0);

	/**@brief The core growing algorithm that will use framents to grow instead of individual atoms
		 @details This feature is still under developement (and thought). It might be implemented or dropped in future
	*/
	void fragmentgrow(ostream& logfile,const System& s,const FragmentSet& frags,int n,int trials=25,int nupper=120,int levellim=LEVELLIMIT,bool restr=false,double autoseed=0,int devid=0)
	{
		int t=0;
    cout << "fragmentgrow(rigid) begins\n";
    logfile << "fragmentgrow(rigid) begins\n";
    bool acc=true;
    Atom* gP=nullptr;
		MolecularFragment* accfrag;
    while(getEffectiveSize()<n)
    {
      acc=true;
      startTrials(devid);
  		if(t>trials)
			{
				cout << "Too many failures: \n"; logfile << "Too many failures: \n";
				break;
			}
  		gP=getRandomBondableAtom();
      if(!gP)
			{
				cout << "Out of valencies\n";
				logfile << "Ended: Out of valencies\n";
				commit();
				break;
			}
      cout << "New begin: "; describeStructure();
			//completeAttemptFragment(ostream& logfile,const System& s,Atom* src,const FragmentSet& fs,int trials,int nupper,int level,int levellim,bool restr,int common_seedlength,int devid)
			bool nextstep=false;
			if(RNG::throwarandompoint(0,1)<autoseed) nextstep=completeAttempt(logfile,s,gP,ATTEMPTTRIALS,nupper,0,levellim,restr,false,false,false,NULLMOL,0,devid);
			else
			{
				accfrag=completeAttemptFragment(logfile,s,gP,frags,ATTEMPTTRIALS,nupper,0,levellim,restr,devid);
				nextstep=(accfrag)?true:false;
				//if(!nextstep) {this->dumpMol("nofragleft.pdb",s.ff); nextstep=completeAttempt(logfile,s,gP,ATTEMPTTRIALS,nupper,0,levellim,restr,false,false,false,NULLMOL,0,devid);} //TODO: Optionally remove this forced DeNovo addition
				if(nextstep) componentFrags.push_back(accfrag);
			}
      if(nextstep)
			{
				cout << "Commiting: "<<gP->toString()<<"\n";
				commit(); t=0;
			}
      else
      {
        cout << "Rollback\n";
        rollback(s);
				#ifdef REALTIME
        	this->dumpMol("rolledback.pdb",s.ff);
					this->dumpMol("realtime.pdb",s.ff);
					usleep(REALTIME*2);
				#endif
        t++;
        cout << "Rolled back: ";
        logfile << "Rolled back: ";
        this->describeStructure(logfile);
      }
    }
		if(this->getEffectiveSize()<n-2) return; //Might not be needed for fragment-based case
    //this->dumpMol("pre_reduction.pdb");

		// The simple physics-based optimization is below (same as Molecule::samplegrow())
		if(USE_PHYSICS)
		{
			this->dumpMol("preopt_fraggrown.pdb",s.ff);
			physics::optimizeMolecule(this,*s.ff,250);
			cout << "Short Geomtery-optimized molecule!\n";
		}
		try
		{
			this->reduce(s,trials,true,devid);
			if(FORCE_CHARGES) this->uniformChargeAdjustments(this->getNetCharge(*s.ff));
		}
		catch(AtomNotInMoleculeException& ex) {this->dumpMol("pre_reduction.pdb"); std::cout << "Crashed with AtomNotInMoleculeException\n"; throw std::exception();}
		if(USE_PHYSICS)
		{
			//physics::gradientDescentPhysics(this);
			physics::optimizeMolecule(this,*s.ff);
			cout << "Well Geomtery-optimized molecule!\n";
		}
	}
	#endif
	/**@brief Sets the net charge of the Molecule by uniformly modifying the charges of all atoms (usually by a small quantity for large enough molecules)
		@details Some atom-types are not perfect in the charge values (especially when repuroposed in two different cases). In this cases there is a minor deviation from the usual charge (0 for neutral molecule), etc.<br/> This minor deviation can build up over multiple such atom-types placed over the molecule. This method corrects for it.
	*/
	void uniformChargeAdjustments(double tchg=0)
	{
		double chg=0; for(int i=0;i<atoms.size();i++) chg+=atoms[i]->seek_charge();
		if(!QUIET)
		{
			cout << "Net charge on system: "<<chg<<"\n";
			cout << "Setting net charge on system to "<<tchg<<" by uniform charge modification\n";
		}
		chg=(tchg-chg)/atoms.size();
		for(int i=0;i<atoms.size();i++) atoms[i]->setCharge(atoms[i]->seek_charge()+chg);
	}
	/**@brief Get expected charge of a molecule.
		 @details The expected charge is calculated by summing the non-zero formal charges on each atom. The formal charges are defined separately in the forcefield (See algorithm PDF)
	*/
	inline double getNetCharge(const ForceField& ff) const
	{
		double nc=0;
		for(Atom* at : atoms) nc+=at->getFormalCharge();
		return nc;
	}

	/**@brief Symbolically describe the structure of the molecule
		 @details Describe the molecule atom-wise along with the free valecies in brackets.<br/>Atoms that must cyclize have a '*' on the side and atoms that have cyclized have a '$'
	*/
	void describeStructure(std::ostream& outs=std::cout) const
	{
		for(Atom* a : atoms) outs << a->toString()<<((a->mustCycl)?"*":"")<<((a->isCyclized())?"$":"")<< "("<<a->seek_valency()<<")";
		outs << "\t";
		for(Atom* a : freeatoms) outs << a->toString() << "("<<a->seek_valency()<<")";
		outs <<"\n";
	}

	/**@brief Locate the closest path joining the atoms at the two indices without passing through the given list of atoms (Used to find path when cyclization occurs)
		 @param[in] latind: Target atom index
		 @param[in] startatin: Starting atom index
		 @param[in] passed: vector of atoms which must not be in the path.
	*/
	std::vector<int> quickLocate(int latind,int startatin,std::vector<int> passed)
	{
		//std::vector<int> passed=passin;
		passed.push_back(startatin);
		for(const int& j : bonds[startatin])
		{
			if(j==latind) return passed;
			if(contains(passed,j)) continue;
			try {return quickLocate(latind,j,passed);}
			catch(NoPathPossibleException& ex) {continue;}
		}
		throw NoPathPossibleException();
	}

	//Output from Molecule to file
	/**@brief Write (dump) a molecule to a file (in PDB format)
		 @details This functions writes the molecule (with atom-type names and coordinates) to a PDB file. A ForceField is also takaen as input (optionally) to detect elements separately from atom-type names (eg. CLAR is chlorine and not carbon)
	*/
	void dumpMol(const std::string& fn,const ForceField* ff=nullptr) const
	{
		cout << "Dumping to: "<<fn<<"\n";
		ofstream of; of.open(fn,ios::out);
		dumpMol(of,ff);
		of.close();
	}
	/**@brief Same as dumpMol(const std::string&,const ForceField*). Uses an output stream instead of filename. */
	void dumpMol(ofstream& f,const ForceField* ff=nullptr) const; //assume open

	/**@brief Dump the molecule as a TXT file (See algorithm documentation PDF for details about this format)
		 @details This function dumps in the old TXT format. This function is kept as backward support.
	*/
	void dumpMolAsTxt(ostream& f,bool cpt=false) const
	{
		f<<atoms.size()<<"\n"; for(int i=0;i<atoms.size();i++) f << *(atoms[i]) << "\n";
		if(cpt)
		{
			for(int i=0;i<atoms.size();i++) {for(const int bmi : bonds[i]) f <<" "<< bmi; f << "\n";}
			f << this->myBE <<" "<<this->myHBE<<"\n";
		}
	}
	/** Same as dumpMolAsTxt(const std::string&). Uses an output stream instead of filename. */
	inline void dumpMolAsTxt(const std::string& f,bool cpt=false) const {ofstream of; of.open(f,ios::out); this->dumpMolAsTxt(of,cpt); of.close();}
	inline void dumpMolAsMol2(const std::string& f,ForceField* ff=nullptr,bool no_orders=false,bool calculate_orders=false,const std::string& name="DNV")
	{
		ofstream of; of.open(f);
		dumpMolAsMol2(of,ff,no_orders,calculate_orders,name);
	}
	/**@brief Write the molecule data in GRO format
		 @details GRO format is the standard GROMACS format.
	*/
	void dumpMolAsGRO(std::ostream& os,int stI=1,bool head=true,const std::string& header="DNV GRO molecule output") const
	{
		if(head)
		{
			os << header << "\n";
			os <<" "<< this->getSize() << "\n";
		}
		char rnum[10],x[10],y[10],z[10],atn[10],ind[10];
		char rnam[10]; sprintf(rnam,"%4s",atoms[0]->getResidue().c_str());
		for(Atom* a : atoms)
		{
			sprintf(rnum,"%d",a->getResidueNumber());
			sprintf(atn,"%+4s",a->toString().c_str());
			sprintf(ind,"%d",stI++);
			sprintf(x,"%.3f",a->seek_x());
			sprintf(y,"%.3f",a->seek_y());
			sprintf(z,"%.3f",a->seek_z());
			char line[100];
			sprintf(line,"%+5s%+5s%+5s%+5s %+7s %+7s %+7s\n",rnum,rnam,atn,ind,x,y,z);
			os << line;
		}
	}
	/**@brief Dump molecule in mol2 format to a stream
		 @details Provides a function write the molecule in mol2 format.
		 @param[in] os: The output stream to write to.
		 @param[in] ff: The ForceField to use for guessing bond-types (and orders)
		 @param[in] no_orders: Ignores orders and bond-type data in mol2 bond records.
		 @param[in] calculate_orders: Calculate order before writing if not already calculated. Must provide ff
		 @param[in] name: Name to give the molecule (Default: DNV)
	*/
	void dumpMolAsMol2(std::ostream& os,ForceField* ff=nullptr,bool no_orders=false,bool calculate_orders=false,const std::string& name="DNV")
	{
		if(no_orders || !ff) calculate_orders=false;
		std::string globalrecs="@<TRIPOS>MOLECULE\n"+name+"\n",atomrecs="@<TRIPOS>ATOM\n",bondrecs="@<TRIPOS>BOND\n";
		int bondcount=0,conK=0;
		const Category* aromcat=nullptr;
		const Category* halocat=nullptr;
		if(ff)
		{
			aromcat=ff->hasCategory("aromatic")?&ff->getCategory("aromatic"):nullptr;
			halocat=ff->hasCategory("halogen")?&ff->getCategory("halogen"):nullptr;
		}
		if(!orders.size() && calculate_orders) this->calculateBondOrders3(*ff);
		for(int i=0;i<atoms.size();i++)
		{
			char atstr[100];
			char fch=atoms[i]->toString()[0];
			std::string symb=(halocat && halocat->contains(atoms[i]))?((fch=='C'?"Cl":(fch=='B'?"Br":std::string(1,fch)))):std::string(1,fch);
			if(aromcat && aromcat->contains(atoms[i]))
			{
				if(fch=='O') symb+=".2";
				else symb+=".ar";
			}
			else
			{
				if(atoms[i]->getHybridization()) symb+="."+to_string(atoms[i]->getHybridization());
			}
			char x[10],y[10],z[10];
			sprintf(x,"%4.4f",atoms[i]->seek_x()*10);
			sprintf(y,"%4.4f",atoms[i]->seek_y()*10);
			sprintf(z,"%4.4f",atoms[i]->seek_z()*10);
			sprintf(atstr,"   %4d %4s     %9s %9s %9s %5s   1  DNV        %1.4f\n",i+1,atoms[i]->toString().c_str(),x,y,z,symb.c_str(),atoms[i]->seek_charge());
			atomrecs+=atstr;
			bondcount+=bonds[i].size();
			for(int j=0;j<bonds[i].size();j++)
			{
				if(bonds[i][j]<i) continue;
				char bstr[60];
				std::string btype="1";
				if(!no_orders)
				{
					if(aromcat && aromcat->contains(atoms[i]) && aromcat->contains(atoms[bonds[i][j]])) btype="ar";
					else {if(orders.size() && orders[i][j]!=1) btype=to_string(orders[i][j]);}
				}
				sprintf(bstr,"  %4d  %4d  %4d  %4s\n",++conK,i+1,bonds[i][j]+1,btype.c_str());
				bondrecs+=bstr;
			}
		}
		globalrecs+=to_string(atoms.size())+" "+to_string(bondcount/2)+" 0 0 0\n";
		if(binding_energy()) os << "#REMARK "<<binding_energy()<<"\t"<<hotspot_binding_energy()<<"\t"<<pE<<"\t"<<errBE<<"\t"<<errHBE<<"(error not impletemented yet -- Last two columns)\n";
		for(const auto& ring : cycles)
		{
			os << "#REMARK CYCLE "<<ring.size()<<" with atoms: ";
			for(Atom* a : ring) os << indexOf(a) <<" ";
			os<< "\n";
		}
		os << globalrecs<<atomrecs<<bondrecs<<"\n";
	}

	/**@brief Get the net solvation energy estimate for this molecule (first-order linear approximation)*/
	float getSolvationEnergyEstimate()
	{
		float solv=0;
		for(Atom* a : atoms) solv+=a->getSolvationCorrection();
		return solv;
	}

	/**@brief Convert adjacency lists into list of bonded pairs. You can choose to repeat each pair twice (Default: No repetition)*/
	std::vector<std::pair<int,int>> generateBondedPairs(bool repeat=false) const
	{
		std::vector<std::pair<int,int>> ret;
		for(int i=0;i<bonds.size();i++)
		{
			for(int j=0;j<bonds[i].size();j++)
			{
				if(!repeat && bonds[i][j]>i) continue;
				ret.push_back(make_pair(i,bonds[i][j]));
			}
		}
		return ret;
	}


  //Geometry
	/**@brief Get the center of geometry of the molecule
		 @details Calculates the center of geometry for a molecule. The position is returned as a Vector (Eigen::Vector3d)
		 @param[in] inch: Include hydrogens? If set to true, the hydrogren atom positions are also considered for calculation of center of geomtery. Default: false
	*/
  Eigen::Vector3d getCentreOfGeometry(bool inch=false) const
  {
    Eigen::Vector3d COG(0,0,0);
    for(Atom* a : getAtoms())
    {
      if(!inch && a->isHydrogen()) continue;
      COG+=a->getPosition();
    }
    if(inch) COG/=getSize();
    else COG/=getEffectiveSize();
    return COG;
  }
	/**@brief Get the center-of-mass for this molecule. Important when handling the system for mechanics*/
	Eigen::Vector3d getCentreOfMass(bool inch=true) const
	{
		Eigen::Vector3d COM(0,0,0);
		double mass=0;
    for(Atom* a : getAtoms())
    {
      if(!inch && a->isHydrogen()) continue;
      COM+=a->getPosition()*a->seek_mass();
			mass+=a->seek_mass();
    }
    COM/=mass;
    return COM;
	}
	#ifdef MACROMODE
	inline Eigen::Vector3d getStoredCentreOfMass(bool inch=true)
	{
		if(comtype!=inch || !hasCOMstored) {storedCOM=getCentreOfMass(inch); hasCOMstored=true; comtype=inch;}
		return storedCOM;
	}
	inline Eigen::Vector3d getStoredCentreOfGeometry(bool inch=true)
	{
		if(cogtype!=inch || !hasCOGstored) {storedCOG=getCentreOfGeometry(inch); hasCOGstored=true; cogtype=inch;}
		return storedCOG;
	}
	inline void storeNormalVector(const Eigen::Vector3d& v) {storedNormal=v;}
	inline const Eigen::Vector3d& getStoredNormalVector() const {return storedNormal;}
	#endif
	/**@brief Calculate the moment of inertia for this molecule (used for simulating rigid-body dynamics)*/
	double getMomentOfInertia(bool inch=true) const
	{
		Eigen::Vector3d COM=getCentreOfMass(inch);
		double MOI=0;
		for(Atom* a : getAtoms())
    {
      if(!inch && a->isHydrogen()) continue;
			Eigen::Vector3d R=(a->getPosition()-COM); double r=R.dot(R);
      MOI+=a->seek_mass()*r*r;
    }
		return MOI;
	}
	/**@brief Get a containing "box" co-ordinates which ensures that the molecule is enclosed by it
		 @details This functions will give the co-ordinates (x<sub>min</sub>,y<sub>min</sub>,z<sub>min</sub>) and (x<sub>max</sub>,y<sub>max</sub>,z<sub>max</sub>).<br/> Using these coordinates, a bounding box can be constructed
		 @return Pair of Vectors (std::pair<Eigen::Vector3d,Eigen::Vector3d>), one for each co-ordinate
	*/
  std::pair<Eigen::Vector3d,Eigen::Vector3d> getContainer() const
  {
    Eigen::Vector3d minC,maxC,temp;
    if(getSize()<=0) return make_pair(Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,0,0));
    minC=atoms[0]->getPosition(); maxC=atoms[0]->getPosition();
    for(int i=1;i<atoms.size();i++)
    {
      temp=atoms[i]->getPosition();
      minC(0)=min(temp(0),minC(0));
      minC(1)=min(temp(1),minC(1));
      minC(2)=min(temp(2),minC(2));

      maxC(0)=max(temp(0),maxC(0));
      maxC(1)=max(temp(1),maxC(1));
      maxC(2)=max(temp(2),maxC(2));
    }
    return make_pair(minC,maxC);
  }
	/**@brief Get the sum of squares distances of each atom (including hydrogen atoms) of this molecule to a given plane
		 @details The plane is taken as input in form of a vector (a,b,c,d) which corresponds to the plane "<b>a</b>x +<b>b</b>y +<b>c</b>z = <b>d</b>".<br/> It is a static function which requires the global variable, "Molecule* currmol" to be set to point to the molecule for which this calculation is to be done.
	*/
  static double squareDistance(std::vector<double> plane)
  {
    double d=0,n=plane[0]*plane[0]+plane[1]*plane[1]+plane[2]*plane[2];
    for(Atom* a : currMol->atoms) d+=sqr(plane[0]*a->seek_x()+plane[1]*a->seek_y()+plane[2]*a->seek_z()-plane[3]);
    return d/n;
  }
	/**@brief Get the bast plane of fit to this molecule - experimental
		 @details This function uses Gradient Descent Fitting to fit a plane equation to this molecule. It assigns the global "currmol" to "this" molecule.<br/>See also: squareDistance(). It calculates the square distance from the plane (which is the target function to be minimized by gradient descent)
		 @param [in] err: The tolerance limit for convergence of gradient descent
		 @return The coeffients of a plane in vector form (a,b,c,d) which corresponds to the plane "<b>a</b>x +<b>b</b>y +<b>c</b>z = <b>d</b>". (std::vector<double>)
	*/
  std::vector<double> bestPlaneOfFit(double err=1e-5) const {currMol=this; return quickmath::gradientDescentFit(squareDistance,4,err);}
	/**@brief Translate a molecule along a given vector. Length of translation depends on vector length*/
	inline void translateMolecule(const Eigen::Vector3d& v) {for(Atom* a : atoms) a->setPosition(a->getPosition()+v);}
	/**@brief Rotate the molecule about a given vector by a given angle (in radians). The angle is given by the norm of the vector
		 @details The rotation is by default around the COM. However, you can choose another center, which will remain fixed during rotation.<br/>
		 If the COM (i.e. centre for the rotation) is chosen to be an atom's position, that atom remains in the same position (useful for rotating around pivotal atoms).
	*/
	void rotateMolecule(const Eigen::Vector3d& n,Eigen::Vector3d com=Eigen::Vector3d(0,0,0))
	{
		if(com(0)==0 && com(1)==0 && com(2)==0) com=this->getCentreOfMass();
		auto rmat=quickgeom::getRotationMatrix(n/n.norm(),n.norm());
		for(Atom* a : atoms)
		{
			Eigen::Vector3d r=a->getPosition()-com;
			r=rmat*r;
			a->setPosition(com+r);
		}
	}
	/**@brief Rotate the dihedral to set it to a particular value
		 @details Provided that the dihedral you have chosen is a bridge bond, the dihedral is rotated (via the smaller side of the dihedral) to match the given value.<br/>
		 The algorithm malfunctions if the chosen bond is part  of a cycle (i.e. not free to rotate). See Molecule::getAllBridgeHeads() to find bridge-head atoms to use for finding freely rotable dihedrals.
	*/
	bool rotateDihedralTo(Atom* a1,Atom* a2,Atom* a3,Atom* a4,double dih)
	{
		int i1=-1,i2=-1,i3=-1,i4=-1;
		for(int i=0;i<atoms.size();i++)
		{
			if(atoms[i]==a1) i1=i;
			else if(atoms[i]==a2) i2=i;
			else if(atoms[i]==a3) i3=i;
			else if(atoms[i]==a4) i4=i;
		}
		if(i1>=0 && i2>=0 && i3>=0 && i4>=0) return rotateDihedralTo(i1,i2,i3,i4,dih);
		else {cout << "In rotateDihedralTo(): Atom chosen is not in given Molecule\n"; cout << i1<<","<<i2<<","<<i3<<","<<i4<<"\n"; throw AtomNotInMoleculeException();}
	}
	/**@brief rotateDihedralTo() that uses atom indices instead of Atom pointers. More efficient if you know the indices already*/
	bool rotateDihedralTo(int i1,int i2,int i3,int i4,double dih)
	{
		std::pair<std::vector<Atom*>,std::vector<Atom*>> sides=splitByBridge(i2,i3);
		double adih=chemtools::getDihedral(atoms[i1],atoms[i2],atoms[i3],atoms[i4]);
		Eigen::Vector3d orient=atoms[i3]->getPosition()-atoms[i2]->getPosition(); orient*=(adih-dih)/orient.norm();
		Molecule* tempm=nullptr; Atom* refa;
		bool ch=get<0>(sides).size()<get<1>(sides).size();
		if(ch) {tempm=new Molecule(get<0>(sides)); refa=atoms[i2];}
		else {tempm=new Molecule(get<1>(sides)); refa=atoms[i3];}
		if(!ch) tempm->rotateMolecule(orient,refa->getPosition());
		else tempm->rotateMolecule(-orient,refa->getPosition());
		if(ch) {for(int i=0;i<tempm->getSize();i++) (get<0>(sides))[i]->setPosition(tempm->getAtom(i)->getPosition());}
		else {for(int i=0;i<tempm->getSize();i++) (get<1>(sides))[i]->setPosition(tempm->getAtom(i)->getPosition());}
		delete tempm;
		return true;
	}

	/**@brief Calculate the total internal energy of a molecule
		 @details Calculates the sum all all forms of internal potential for a molecule. It includes angle, dihedral, and non-bonding interaction energies.<br/>This function still requires proper testing.
	*/
	double calculateInternalEnergy(const ForceField& ff) const
	{
		double inte=0;
		double ange=0,dihe=0,ofe=0,lje=0;
		std::vector<Atom*> visited;
		for(int i=0;i<atoms.size();i++)
		{
			if(atoms[i]->isHydrogen()) continue; //Comment to calculate H-H interactios as well.
			visited.push_back(atoms[i]);
			std::vector<Atom*> bndats=this->getBondedAtoms(atoms[i]);
			std::vector<Atom*> nearat=bndats;
			for(Atom* ca : bndats)
			{
				if(ca->isHydrogen()) continue;
				std::vector<Atom*> sbndats=this->getBondedAtoms(ca);
				for(Atom* sca : sbndats)
				{
					if(sca->isHydrogen() || sca==atoms[i]) continue;
					nearat.push_back(sca);
					std::vector<Atom*> tbndats=this->getBondedAtoms(ca);
					for(Atom* tca : tbndats)
					{
						if(tca==atoms[i]) continue;
						if(tca->isHydrogen()) continue; //Comment to calculate 1-4s for Hydrogen as well
						ofe+=chemtools::getNonbondingPotential(tca,atoms[i])/2;
						nearat.push_back(tca);
					}
				}
			}
			for(Atom* fa : atoms)
			{
				if(fa->isHydrogen()) continue; //Comment to calculate LJs for hydrogen as well
				if(fa==atoms[i] || contains(nearat,fa)) continue;
				double ae=chemtools::getNonbondingPotential(fa,atoms[i]);
				cout << atoms[i]->toString()<<" "<<atoms[i]->getPosition().transpose() << " "<<fa->toString()<<" "<<fa->getPosition().transpose() << "\t"<<ae<<"\n";
				lje+=ae;
			}
		}
		ange/=2; dihe/=2; lje/=2; ofe/=2;
		cout << "Angle energy: "<<ange<<"\n";
		cout << "Dihedral energy: "<<dihe<<"\n";
		cout << "LJ Energy: "<< lje<<"\n";
		cout << "1-4 Energy: "<<ofe<<"\n";
		return ange+dihe+lje+ofe;
	}
	/**@brief Calculate the forces acting on a target molecule from this molecule*/
	std::pair<std::vector<Eigen::Vector3d>,std::vector<Eigen::Vector3d>> calculateActingForces(Molecule* lig) const
	{
		std::vector<Eigen::Vector3d> forces,poses;
		Eigen::Vector3d com=lig->getCentreOfMass();
		for(Atom* a : lig->getAtoms())
		{
			Eigen::Vector3d frc(0,0,0);
			for(Atom* p : this->getAtoms())
			{
				Eigen::Vector3d r=p->bondVectorFrom(a); double n=r.norm();
				frc+=(chemtools::getLJForce(p,a,n)/n)*r+(chemtools::getElectrostaticForce(p,a,n)/n)*r;
			}
			poses.push_back((a->getPosition()-com).cross(frc));
			forces.push_back(frc);
		}
		return make_pair(forces,poses);
	}

	inline bool inSameCycle(Atom* a1,Atom* a2) const
	{
		bool f1=false,f2=false;
		for(auto cd : cycles)
		{
			f1=false;f2=false;
			for(Atom* tg : cd)
			{
				if(tg==a1) {f1=true; if(f2) break;}
				if(tg==a2) {f2=true; if(f1) break;}
			}
			if(f1 && f2) return true;
		}
		if(bonds[this->indexOf(a2)].size()<2) return false;
		std::vector<Atom*> catt; catt.push_back(a1);
		int paths=0;
		Atom* foc=a1;
		std::vector<Atom*> nxtfocs;
		while(foc)
		{
			catt.push_back(foc);
			std::vector<Atom*> nnh=this->getBondedAtoms(foc);
			for(Atom* fn : nnh)
			{
				if(fn->getStandardValency()<2 || contains(catt,fn)) continue;
				if(fn==a2) {paths++; if(paths==2) break;}
				else
				{
					nxtfocs.push_back(fn);
				}
			}
			if(nxtfocs.size()) {foc=nxtfocs[nxtfocs.size()-1]; nxtfocs.pop_back();}
			else foc=nullptr;
		}
		return paths>=2;
	}

	/**@brief Find all bridge-head atom index pairs for all possible bridges in this molecule
		 @details Looks for bridge edges in the molecule and produces all pairs of atoms that have bridge-bonds.
	*/
	std::vector<std::pair<int,int>> getAllBridgeHeads()
	{
		std::vector<std::pair<int,int>> ret;
		this->markCyclizedAtoms();
		for(int i=0;i<atoms.size();i++)
		{
			if(atoms[i]->getStandardValency()<2) continue;
			for(int j=0;j<bonds[i].size();j++)
			{
				if(bonds[i][j]>i) continue;
				if(atoms[bonds[i][j]]->getStandardValency()<2) continue;
				if(!inSameCycle(atoms[i],atoms[bonds[i][j]])) ret.push_back(make_pair(i,bonds[i][j]));
			}
		}
		return ret;
	}
	std::pair<std::vector<Atom*>,std::vector<Atom*>> splitByBridge(int a1,int a2)
	{
		if(!contains(bonds[a1],a2)) throw NoSuchBondException();
		std::vector<Atom*> f1,f2;
		std::vector<int> visited(1,a1);
		bool end=false;
		while(!end)
		{
			end=true;
			int sl=visited.size();
			for(int j=0;j<sl;j++)
			{
				int ci=visited[j];
				for(int i=0;i<bonds[ci].size();i++)
				{
					if(bonds[ci][i]==a2 || contains(visited,bonds[ci][i])) continue;
					visited.push_back(bonds[ci][i]);
					end=false;
				}
			}
		}
		for(int i : visited) f1.push_back(atoms[i]);
		end=false;
		visited=std::vector<int>(); visited.push_back(a2);
		while(!end)
		{
			end=true;
			int sl=visited.size();
			for(int j=0;j<sl;j++)
			{
				int ci=visited[j];
				for(int i=0;i<bonds[ci].size();i++)
				{
					if(bonds[ci][i]==a1 || contains(visited,bonds[ci][i])) continue;
					visited.push_back(bonds[ci][i]);
					end=false;
				}
			}
		}
		for(int i : visited) f2.push_back(atoms[i]);
		return make_pair(f1,f2);
		//return make_pair(createSubmolecule(f1),createSubmolecule(f2));
	}

private:
	std::vector<Atom*> parseForSMILES2(int atind,std::string& smiles,std::vector<std::string>& atomnames,const std::vector<int>& hcounts,std::vector<int>& visited,std::vector<Atom*>& atorder,int& ringID,int fromind=-1) const
	{
		visited.push_back(atind);
		atorder.push_back(atoms[atind]);
		smiles+="&"+to_string(atind)+"$";
		int endI=bonds[atind].size()-1;
		for(int i=0;i<bonds[atind].size();i++) {if(atoms[bonds[atind][i]]->isHydrogen()) continue; endI=i;}
		for(int i=0;i<bonds[atind].size();i++)
		{
			int cI=bonds[atind][i];
			if(cI==fromind) continue;
			if(atoms[cI]->isHydrogen()) continue;
			if(contains(visited,cI))
			{
				if(cI>atind) continue;
				//cout << "It is already visited: "<<atoms[cI]->toString()<<"("<<cI<<")-"<<atoms[atind]->toString()<<"("<<i<<")"<<"\n";
				if(ringID>10) {atomnames[cI]+="%"; atomnames[atind]+="%";}
				atomnames[cI]+=to_string(ringID);
				switch(orders[atind][i])
				{
					case 2:
						atomnames[atind]+="=";
						break;
					case 3:
						atomnames[atind]+="#";
						break;
				}
				atomnames[atind]+=to_string(ringID++);
			}
			else
			{
				if(i!=endI) smiles+="(";
				switch(orders[atind][i])
				{
					case 2:
						smiles+="=";
						break;
					case 3:
						smiles+="#";
						break;
				}
				parseForSMILES2(cI,smiles,atomnames,hcounts,visited,atorder,ringID,atind);
				if(i!=endI) smiles+=")";
			}
		}
		return atorder;
	}
public:
	/**@brief This function allows to drop one atom from the molecule. Read well before usage
		 @details Using this function, you can drop any atom which is a terminal atom. This is most commonly a hydrogen atom, and sometimes a halogen, but more uncommonly it might also be carbonyl oxygen.<br/>You can also remove other atoms as long as you are sure they have exactly one connection to them at the time of removal.
		 @param[in] din: Index of atom to be removed
	*/
	void dropSingleAtom(int din)
	{
		std::vector<int> conA=bonds[din];
		for(int& c : conA)
		{
			int i=0;
			for(;i<bonds[c].size();i++) {if(bonds[c][i]==din) break;}
			bonds[c].erase(bonds[c].begin()+i);
			atoms[c]->setValency(atoms[c]->seek_valency()+1);
			if(!contains(freeatoms,atoms[c])) freeatoms.push_back(atoms[c]);
		}
		bonds.erase(bonds.begin()+din);
		myBE-=atoms[din]->getEnergyContribution();
		myHBE-=atoms[din]->getHotspotEnergyContribution();
		if(atoms[din]->isHydrogen()) HC--;
		if(VARIEDCOUNT && VCLIST->contains(atoms[din])) HAC--;
		atoms.erase(atoms.begin()+din);
		for(std::vector<int>& bv : bonds)
		{
			for(int& i : bv) {if(i>=din) i--;}
		}
	}
	/**@brief This function generates the starting fragment for "pruning" (See algorthim PDF).
		 @details The algorithm of pruning requires removal of random "terminal" atoms. Atoms which have a single non-hydrogen neighbour are treated as terminals, and the non-hydrogen neighbour is called the "stem".
		 @param[in] n: It is the prune depth. 'n' terminal atoms are recursively removed before pruning begins. After the first terminal atom is removed, the next atom (which initially was not terminal but is terminal now) can be removed (hence the term 'depth')
	*/
	void makePruneFragment(int n,int devid=0)
	{
		cout << "Begin. n="<<n<<"\n";
		std::vector<std::vector<int>> hydats;
		std::vector<int> terms,termconns;
		bool terminal=true;
		int termholder=-1;
		std::vector<int> hydatis;
		for(int i=0;i<atoms.size();i++)
		{
			if(atoms[i]->isHydrogen()) {hydatis.push_back(i); continue;}
			std::vector<int> hyd;
			terminal=true;
			termholder=-1;
			for(const int& bai : bonds[i])
			{
				if(atoms[bai]->isHydrogen()) hyd.push_back(bai);
				else
				{
					if(termholder==-1) termholder=bai;
					else {terminal=false; break;}
				}
			}
			if(terminal)
			{
				if(termholder==-1) cout << "WARN: holder not found!!!\n";
				terms.push_back(i); hydats.push_back(hyd); termconns.push_back(termholder);
			}
		}
		if(!terms.size())
		{
			if(hydatis.size())
			{
				cout << "NOTE: Prune fragment generation will only remove hydrogen atoms now. No heavy-atom terminals to prune.\n";
				int rhi=randomSelect(hydatis,devid);
				this->dropSingleAtom(rhi);
				n--;
				if(n>0) makePruneFragment(n,devid);
				else this->dumpMol("prunefrag.pdb");
				return;
			}
			else {cout << "WARN: Prune fragment generation failed. No terminals to prune. Stopping.\n"; describeStructure(); return;}
		}
		cout << "Terminals:\n"; for(int i=0;i<terms.size();i++) cout<<terms[i]<<":"<< atoms[terms[i]]->toString()<<"\n";
		int rS=(int)RNG::throwarandompoint(devid,0,terms.size());
		cout << "Chose rS="<<rS<<" i.e. "<<terms[rS]<<":"<<atoms[terms[rS]]->toString()<<"\n";
		std::vector<int> erasure(1,terms[rS]); for(const int& hai : hydats[rS]) erasure.push_back(hai);
		int i=-1;
		for(i=0;i<bonds[termconns[rS]].size();i++)
		{
			if(bonds[termconns[rS]][i]==terms[rS]) break;
		}
		cout << "Using bonds of "<<termconns[rS]<<":"<<atoms[termconns[rS]]->toString()<<"\n";
		if(i>=bonds[termconns[rS]].size()) cout << "Failed to find bond for "<<terms[rS]<<" in: "<<bonds[termconns[rS]]<<"\n";
		cout << terms[rS] << ":"<<atoms[terms[rS]]->toString()<<"\t"<<i<<"\n";
		bonds[termconns[rS]].erase(bonds[termconns[rS]].begin()+i);
		if(!atoms[termconns[rS]]->seek_valency()) freeatoms.push_back(atoms[termconns[rS]]);
		atoms[termconns[rS]]->setValency(atoms[termconns[rS]]->seek_valency()+1);
		HC-=hydats[rS].size();
		if(VARIEDCOUNT && VCLIST->contains(atoms[terms[rS]])) HAC--;
		for(int i=0;i<erasure.size();i++)
		{
			myBE-=atoms[erasure[i]]->getEnergyContribution();
			myHBE-=atoms[erasure[i]]->getHotspotEnergyContribution();
			for(int deli=0;deli<freeatoms.size();deli++) {if(freeatoms[deli]==atoms[erasure[i]]) {freeatoms.erase(freeatoms.begin()+deli);break;}}
			delete atoms[erasure[i]];
			atoms.erase(atoms.begin()+erasure[i]);
			bonds.erase(bonds.begin()+erasure[i]);
			for(int p=i+1;p<erasure.size();p++) {if(erasure[p]>erasure[i]) erasure[p]--;}
			for(std::vector<int>& bv : bonds)
			{
				for(int j=0;j<bv.size();j++) {if(bv[j]>erasure[i]) bv[j]--;}
			}
		}
		describeStructure();
		n--;
		if(n>0) makePruneFragment(n);
		else this->dumpMol("prunefrag.pdb");
	}
	inline bool makeRandomPruneFragment(const ForceField* ff, float pcent,bool keepHs=false,int devid=0) {return get<0>(this->makeRandomPruneFragmentWithMapping(ff,pcent,keepHs,devid));}
	std::pair<bool,std::map<int,int>> makeRandomPruneFragmentWithMapping(const ForceField* ff, float pcent,bool keepHs=false,int devid=0)
	{
		if(pcent>1) pcent/=100; //Assume percentage is given out of 100%
		assert(pcent<1);
		std::map<int,int> mapping;
		std::cout << "XXX - Atoms: "+std::to_string(atoms.size())+" and bonds: "+std::to_string(bonds.size())+"\n";
		int rsize=this->getEffectiveSize()*pcent;
		if(rsize<=1) return std::make_pair(false,mapping);
		std::vector<int> selats;
		int seed=-1;
		while(seed<0)
		{
			seed=RNG::throwarandompoint(devid,0,1)*atoms.size();
			if(atoms[seed]->isHydrogen()) {seed=-1; continue;}
		}
		selats.push_back(seed);
		while(selats.size()<rsize)
		{
			int nxt=seed;
			while(contains(selats,nxt))
			{
				int node=randomSelect(selats);
				nxt=randomSelect(bonds[node]);
				if(atoms[nxt]->isHydrogen()) {nxt=seed; continue;}
			}
			selats.push_back(nxt);
		}
		if(keepHs)
		{
			for(int idx=0;idx<selats.size();idx++)
			{
				int ix=selats[idx];
				for(int hi : bonds[ix])
				{
					if(atoms[hi]->isHydrogen()) selats.push_back(hi);
				}
			}
		}
		for(int i=0;i<selats.size();i++) mapping[i]=selats[i];
		std::vector<Atom*> leftatoms; for(int ix : selats) leftatoms.push_back(atoms[ix]);
		std::vector<std::vector<int>> nbonds;
		for(int i=0;i<selats.size();i++)
		{
			nbonds.push_back(std::vector<int>());
			int ii=selats[i];
			for(int j=0;j<bonds[ii].size();j++)
			{
				int tidx=bonds[ii][j];
				for(int k=0;k<selats.size();k++)
				{
					if(tidx==selats[k]) {nbonds[i].push_back(k); break;}
				}
			}
		}
		atoms=leftatoms;
		bonds=nbonds;
		valencies=std::vector<int>();
		for(int i=0;i<bonds.size();i++)
		{
			int effval=atoms[i]->getStandardValency()-bonds[i].size();
			atoms[i]->setValency(effval);
			valencies.push_back(effval);
		}
		this->HAC=this->getEffectiveSize();
		this->HC=0;
		bool checkcyc=ff->hasCategory("ring");
		for(Atom* at : atoms)
		{
			if(at->isHydrogen()) HC++;
			if(checkcyc && ff->getCategory("ring").contains(at)) at->mustCycl=true;
		}
		this->calculateFreeAtoms();
		if(checkcyc) this->markCyclizedAtoms();
		return std::make_pair(true,mapping);
	}

	/**@brief Make a sub-Molecule from the given molecule using a set of atoms (maintaining any bonds possible)*/
	Molecule* createSubmolecule(const std::vector<Atom*>& ats)
	{
		std::vector<int> relpos; std::vector<std::vector<int>> bdset; for(int i=0;i<ats.size();i++) bdset.push_back(std::vector<int>());
		int backrelpos[atoms.size()+1]; for(int i=0;i<=atoms.size();i++) backrelpos[i]=-1;
		for(int i=0;i<ats.size();i++)
		{
			int aI=this->indexOf(ats[i]);
			relpos.push_back(aI);
			backrelpos[aI]=i;
		}
		//cout << "Relative positions established\n";
		for(int i=0;i<relpos.size();i++)
		{
			for(int j=0;j<bonds[relpos[i]].size();j++)
			{
				if(backrelpos[bonds[relpos[i]][j]]!=-1) bdset[i].push_back(backrelpos[bonds[relpos[i]][j]]);
			}
		}
		std::vector<Atom*> finats; for(Atom* a : ats) finats.push_back(new Atom(*a));
		Molecule* ret=new Molecule(finats,bdset);
		return ret;
	}
	/**@brief Get SMILES format for this molecule - Improved
		 @details This function parses the graph of the molecule to generate a valid smiles string. It requires a ForceField (in order to get valency and hybridization).<br/>For complete molecules (with no free valecies), the "empty" force-field can be used (See algorithm PDF).<br/>It can also be called without a ForceField, in which case it guesses at the valencies using standard parameters.<br/>This function requires the calculateBondOrders() to be called eariler on this Molecule object
		 @param[in] ff: ForceField to use (Note that it's a pointer to the ForceField object)
		 @param[in] daylight: Generate daylight SMILES (represents aromatic atoms with small letters and omits their bond-orders)
	*/
	std::pair<std::string,std::vector<Atom*>> getSMILES2(const ForceField* ff=nullptr,bool daylight=false) const
	{
		if(daylight && (!ff || !ff->hasCategory("aromatic"))) throw NoAromaticCategoryException();
		if(!orders.size()) throw BondOrdersUncalculatedException();
		const Category *aromcat=&(ff->getCategory("aromatic")),*halocat=nullptr; if(ff->hasCategory("halogen")) halocat=&(ff->getCategory("halogen"));
		std::string SMILES="";
		int K=0,ringID=1;
		for(K=0;K<atoms.size();K++) if(!atoms[K]->isHydrogen()) break;
		std::vector<int> hcount,visited;
		for(int i=0;i<atoms.size();i++)
		{
			int hC=0;
			if(atoms[i]->isHydrogen()) {hcount.push_back(0); continue;}
			for(int j=0;j<bonds[i].size();j++) {if(atoms[bonds[i][j]]->isHydrogen()) hC++;}
			hcount.push_back(hC);
		}
		std::vector<std::string> atomnames;
		for(Atom* a : atoms)
		{
			if(halocat && halocat->contains(a))
			{
				if(a->toString()[0]=='C') atomnames.push_back("Cl");
				else if(a->toString()[0]=='B') atomnames.push_back("Br");
				else atomnames.push_back(std::string(1,a->toString()[0]));
			}
			else atomnames.push_back(std::string(1,a->toString()[0]));
			if(a->getFormalCharge())
			{
				int atnind=indexOf(a); if(hcount[atnind]) atomnames[atomnames.size()-1]+="H"+((hcount[atnind]==1)?"":to_string(hcount[atnind]));
				if(a->getFormalCharge()>0) atomnames[atomnames.size()-1]="["+atomnames[atomnames.size()-1]+"+"+to_string((int)a->getFormalCharge())+"]";
				else atomnames[atomnames.size()-1]="["+atomnames[atomnames.size()-1]+to_string((int)a->getFormalCharge())+"]";
			}
		}

		std::vector<Atom*> neworder;
		parseForSMILES2(K,SMILES,atomnames,hcount,visited,neworder,ringID);
		if(!QUIET) cout << "Raw SMILES: "<<SMILES<<"\n";

		//Replacing atom references by actual node data
		for(int i=atoms.size()-1;i>=0;i--)
		{
			if(!atoms[i]->isHydrogen()) cout << atoms[i]->toString() << "\t"<<atomnames[i]<<"\n";
			else continue;
			std::string oldsmiles="";
			while(oldsmiles!=SMILES)
			{
				oldsmiles=SMILES;
				stringfx::replace(SMILES,"&"+to_string(i)+"$",atomnames[i]);
			}
		}

		//Return the result
		if(!QUIET) cout << "New order: \n";
		for(Atom* a : neworder) cout << a->toString() <<"("<<this->indexOf(a)<<")"<< "\n";
		return make_pair(SMILES,neworder);
	}

	//Backup - Rewriting of SMILES generation (still works, mostly)
	/*std::string produceSMILES2(int cai,std::vector<Atom*>& order,std::vector<int>& visited,std::vector<std::string>& names,int ringID=0,const std::vector<int>& branches=std::vector<int>(),int lastI=-1) const
	{
		Atom* ca=atoms[cai];
		order.push_back(ca); visited.push_back(cai);
		std::string smiles=names[cai];
		cout << cai << " named "<<ca->toString()<<"\n";
		int br=0;
		for(int ii=0;ii<bonds[cai].size();ii++)
		{
			if(bonds[cai][ii]==lastI) continue;
			if(atoms[bonds[cai][ii]]->isHydrogen()) continue;
			br++;
			if(contains(visited,bonds[cai][ii]))
			{
				if(cai<bonds[cai][ii]) continue; //Avoid double-labelling
				std::string add=to_string(++ringID);
				if(ringID>=10) add="%"+add;
				names[cai]+=add;
				names[bonds[cai][ii]]+=add;
				continue;
			}
			smiles=smiles+"(";
			switch(orders[cai][ii])
			{
				case 2:
					smiles=smiles+"=";
					break;
				case 3:
					smiles=smiles+"#";
					break;
			}
			smiles=smiles+produceSMILES2(bonds[cai][ii],order,visited,names,ringID,branches,cai);
			smiles=smiles+")";
		}
		return smiles;
	}

	std::pair<std::string,std::vector<Atom*>> getSMILES2(const ForceField* ff=nullptr,bool daylight=false) const
	{
		//Assume bond-orders have been assigned
		const Category *aromcat=&(ff->getCategory("aromatic")),*halocat=nullptr; if(ff->hasCategory("halogen")) halocat=&(ff->getCategory("halogen"));
		std::vector<std::string> atomnames;
		std::vector<int> hcount,visited,branches;
		for(int i=0;i<atoms.size();i++)
		{
			int hC=0;
			if(atoms[i]->isHydrogen()) {hcount.push_back(0); continue;}
			for(int j=0;j<bonds[i].size();j++) {if(atoms[bonds[i][j]]->isHydrogen()) hC++;}
			hcount.push_back(hC);
		}
		for(Atom* a : atoms)
		{
			if(halocat && halocat->contains(a))
			{
				if(a->toString()[0]=='C') atomnames.push_back("Cl");
				else if(a->toString()[0]=='B') atomnames.push_back("Br");
				else atomnames.push_back(std::string(1,a->toString()[0]));
			}
			else if(daylight && aromcat && aromcat->contains(a)) atomnames.push_back(std::string(1,a->toString()[0]+32)); //Lower case
			else atomnames.push_back(std::string(1,a->toString()[0]));
			if(a->getFormalCharge())
			{
				int atnind=indexOf(a); if(hcount[atnind]) atomnames[atomnames.size()-1]+="H"+((hcount[atnind]==1)?"":to_string(hcount[atnind]));
				if(a->getFormalCharge()>0) atomnames[atomnames.size()-1]="["+atomnames[atomnames.size()-1]+"+"+to_string((int)a->getFormalCharge())+"]";
				else atomnames[atomnames.size()-1]="["+atomnames[atomnames.size()-1]+to_string((int)a->getFormalCharge())+"]";
			}
		}
		for(Atom* a : atoms) {if(!a->isHydrogen()) branches.push_back(0);}
		int K=0;
		Atom* stat=atoms[K++]; while(stat->isHydrogen() && K<atoms.size()) stat=atoms[K++];
		visited.push_back(K-1);
		cout << "Starting at K="<<K<<"\n";
		std::vector<Atom*> dumpatorder;
		this->produceSMILES2(K-1,dumpatorder,visited,atomnames,0,branches);
		std::cout << "Round 1 complete\n";
		visited=std::vector<int>();
		std::vector<Atom*> atorder;
		std::string res=this->produceSMILES2(K-1,atorder,visited,atomnames,0,branches);
		return make_pair(res,atorder);
	}*/


	/**@brief Get Atom ID (used during "hashing") -experimental*/
	inline int getAID(Atom* a) const {return getAID(indexOf(a));}
	/**@brief Get Atom ID (used during "hashing") -experimental*/
	inline int getAID(int atno) const {return atoms[atno]->seek_mass()+bonds[atno].size()*3;}
	inline std::vector<int>& hashatomcomplete(Atom* at) const {return hashatomcomplete(indexOf(at));}
	/**@brief Get a "hash vector" for this atom -experimental*/
	std::vector<int>& hashatomcomplete(int ind) const
	{
		//cout << "Hashatom complete called: "<<ind<<"\t"<<atoms[ind]->toString()<<"\n";
		std::vector<int> myhash(1,getAID(ind));
		/*if(bonds[ind].size()<=1)
		{
			if(!atoms[bonds[ind][0]]->hashvector.size()) hashatomcomplete(bonds[ind][0]);
			for(int& i : atoms[bonds[ind][0]]->hashvector) myhash.push_back(i);
			return (atoms[ind]->hashvector=myhash);
		}*/
		std::vector<int> visited(1,ind);
		std::vector<int> depcue(1,ind);
		while(visited.size()<atoms.size())
		{
			std::vector<int> nvis;
			for(const int& ai : depcue)
			{
				for(const int& nai : bonds[ai])
				{
					if(contains(visited,nai)) continue;
					if(bonds[nai].size()<=1)
					{
						visited.push_back(nai);
						break;
					}
					nvis.push_back(nai);
				}
			}
			int dN=0;
			for(const int& ni : nvis) {dN+=getAID(ni);	visited.push_back(ni);}
			myhash.push_back(dN);
			depcue=nvis;
		}
		return (atoms[ind]->hashvector=myhash);
	}
	/**@brief Numeric first-order hash. \deprecated Use hashatomcomplete() now*/
	long int hashatom(Atom* a)
	{
		std::vector<Atom*> bnd=this->getBondedAtoms(a);
		long int myhash=(long int)(a->toString()[0])*10000;
		for(Atom* n : bnd)
		{
			if(bnd.size()==1)
			{
				if(!n->hashvalue) hashatom(n);
				myhash+=n->hashvalue;
			}
			myhash+=(int)(n->toString()[0]*10)*(this->getBondedAtoms(n).size());
			std::vector<Atom*> sobnd=this->getBondedAtoms(n);
			for(Atom* n2 : sobnd)
			{
				int nb=bonds[indexOf(n2)].size();
				if(nb==1) continue;
				if(n2==a) continue;
				myhash+=(long int)(n2->toString()[0])*nb;
			}
		}
		a->hashvalue=myhash;
		return myhash;
	}
	/**@brief Numeric first-order hash. \deprecated Use getCompleteHash() now*/
	long int getFirstOrderHash()
	{
		if(!atoms.size()) return 0;
		long int hash=0;
		for(int i=0;i<atoms.size();i++) {hash+=hashatom(atoms[i]);}
		return hash;
	}
	/**@brief Get a complete hash vector for a molecule. - experimental
		 @details The hash vector generated by this function can uniquely identify each atom of the molecule and the molecule as a whole.<br/>As of now, this hash cannot be used as a fingerprint for similarity matching. However, it can be used to tell if two molecule graphs are exactly same or not
		 @return The hash vector for the entire molecule. If two molecules are same, they will have the same hash vector.
	*/
	std::vector<int> getCompleteHash(bool heavy=false) const
	{
		std::vector<int> fhash;
		for(int i=0;i<atoms.size();i++)
		{
			cout << atoms[i]->toString() << "\n";
			if(!atoms[i]->hashvector.size()) hashatomcomplete(i);
		}
		for(Atom* a : atoms)
		{
			if(heavy && a->isHydrogen()) continue;
			for(int i=0;i<a->hashvector.size();i++)
			{
				if(fhash.size()>i) fhash[i]+=a->hashvector[i];
				else fhash.push_back(a->hashvector[i]);
			}
		}
		return fhash;
	}
	/**@brief Try to match the similarity of two molecules and score the similarity in (0,1) - under development*/
	double getGraphDistance(const Molecule* m2,double base=0.8) const
	{
		std::vector<int> hash1=this->getCompleteHash(true),hash2=m2->getCompleteHash(true);
		int hA=0;
		double baseline=0,avcon=0;
		if(getEffectiveSize()<m2->getEffectiveSize())
		{
			for(int i=0;i<atoms.size();i++)
			{
				if(atoms[i]->isHydrogen()) continue;
				hA++;
				baseline+=atoms[i]->seek_mass();
				avcon+=bonds[i].size();
			}
			baseline/=hA;
			avcon/=hA;
		}
		else
		{
			for(int i=0;i<m2->atoms.size();i++)
			{
				if(m2->atoms[i]->isHydrogen()) continue;
				hA++;
				baseline+=m2->atoms[i]->seek_mass();
				avcon+=m2->bonds[i].size();
			}
			baseline/=hA;
			avcon/=hA;
		}
		cout << "Baseline: "<<(baseline+12)<<"\n";
		cout << "Avcon: "<<avcon<<"\n";
		baseline=(baseline+12)*avcon/2;
		return normalizedVectorDifference(hash1,hash2,this->getEffectiveSize(),m2->getEffectiveSize(),base);
	}
	/**@brief Match to vectors and calculate some kind of a "distance" to score the similarity in (0,1) - under development*/
	static double normalizedVectorDifference(const std::vector<int>& hash1,const std::vector<int>& hash2,int s1,int s2,double base=0.8)
	{
		//std::vector<int> diff;
		double eN=0,nN=(pow(1/base,min(hash1.size(),hash2.size())-1)-1)/(1/base-1);
		for(int i=0;i<min(hash1.size(),hash2.size());i++)
		{
			eN+=(1-(abs(hash1[i]-hash2[i])/max(hash1[i],hash2[i])))*pow(base,-i);
			//if(i!=0) nN+=min(s1,s2)*(baseline)*pow(base,-i); else nN+=min(s1,s2)*(baseline)/4;/* nN+=max(hash1[i],hash2[i])*pow(base,-i);*/
		}
		/*if(hash1.size()>hash2.size()) for(int i=diff.size();i<hash1.size();i++) {diff.push_back(hash1[i]); nN+=hash1[i]*pow(base,-i);}
		else for(int i=diff.size();i<hash2.size();i++) {diff.push_back(hash2[i]); nN+=hash2[i]*pow(base,-i);}*/
		//for(int i=0;i<diff.size();i++) eN+=diff[i]*pow(base,-i);
		cout <<min(hash1.size(),hash2.size())<<":\t"<<s1<<","<<s2<<"\t"<<eN<<","<<nN<<"\n";
		return eN/nN;
	}

	/**@brief Get the "pose" (3D atom-position data) of the Molecule object*/
	inline Pose extractPose()
	{
		std::vector<Eigen::Vector3d> ret;
		for(Atom* a : atoms) ret.push_back(a->getPosition());
		return ret;
	}
	inline void assignPose(const Pose& p)
	{
		for(int i=0;i<p.size();i++) atoms[i]->setPosition(p[i]);
	}
	std::vector<Pose> getAllDihedralPosesFrom(const std::vector<std::vector<Atom*>>& dihats,const ForceField& ff,int sh=0)
	{
		if(!dihats.size()) return std::vector<Pose>(1,this->extractPose());
		std::vector<Pose> ret;
		std::vector<Atom*> curdih=dihats[sh];
		std::vector<double> pars=ff.getDihedralParametersFailsafe(curdih[0],curdih[1],curdih[2],curdih[3]);
		std::vector<double> muldih=ForceField::expandMultiplicity(pars[0],pars[2]);
		cout <<dihats.size()<<":\t"<< sh << " - "<<muldih.size()<<"\n";
		for(int i=0;i<muldih.size();i++)
		{
			this->rotateDihedralTo(curdih[0],curdih[1],curdih[2],curdih[3],muldih[i]);
			if(sh<dihats.size()-1)
			{
				std::vector<Pose> np=getAllDihedralPosesFrom(dihats,ff,sh+1);
				for(Pose p : np) ret.push_back(p);
			}
			else ret.push_back(this->extractPose());
		}
		cout << ret.size()<< " poses finally\n";
		return ret;
	}
};

class MolecularFragment : public Molecule
{
protected:
	std::vector<std::vector<int>> supports;
	std::vector<int> hands;
	Eigen::Vector3d myCOG=Eigen::Vector3d(0,0,0);
	std::vector<Pose> rotatedPoses;
	std::string fragmentlabel="";
public:
	MolecularFragment(std::istream& f,bool eval,bool cpt,const ForceField* ff=nullptr,const std::vector<int>& handles=std::vector<int>(),int prec=10,bool dihman=true) : Molecule(f,eval,cpt,ff)
	{
		hands=handles;
		if(!handles.size()) {for(Atom* a : freeatoms) hands.push_back(indexOf(a));}
		myCOG=this->getCentreOfGeometry();
		calcSupports();
		calcRotatedPoses(*ff,prec,dihman);
	}
	MolecularFragment(const std::string& f,const std::string& fmt,const ForceField& ff,bool eval,const std::vector<int>& handles=std::vector<int>(),int prec=10,bool dihman=true) : Molecule(f,ff,fmt,eval)
	{
		hands=handles;
		if(!handles.size()) {for(Atom* a : freeatoms) hands.push_back(indexOf(a));}
		myCOG=this->getCentreOfGeometry();
		calcSupports();
		calcRotatedPoses(ff,prec,dihman);
	}
private:
	void calcSupports()
	{
		supports=std::vector<std::vector<int>>();
		for(int i=0;i<hands.size();i++)
		{
			std::vector<int> sv;
			for(int& a1i : bonds[hands[i]])
			{
				if(bonds[a1i].size()>1)
				{
					sv.push_back(a1i);
					for(int& a2i : bonds[a1i])
					{
						if(a2i!=hands[i])
						{
							sv.push_back(a2i);
							break;
						}
					}
					break;
				}
			}
			supports.push_back(sv);
		}
	}
	void calcRotatedPoses(const ForceField& ff,int prec=10,bool dihman=true)
	{
		prec++;
		std::vector<std::pair<int,int>> brids=this->getAllBridgeHeads();
		std::vector<std::vector<Atom*>> dihseq;
		if(dihman)
		{
			for(auto p : brids)
			{
				Atom *a2=atoms[get<0>(p)],*a3=atoms[get<1>(p)];
				Atom *a1=nullptr,*a4=nullptr;
				std::vector<Atom*> dihats;
				for(Atom* ax : this->getBondedAtoms(a2))
				{
					if(ax->isHydrogen() || ax==a3) continue;
					else {a1=ax; break;}
				}
				for(Atom* ax : this->getBondedAtoms(a3))
				{
					if(ax->isHydrogen() || ax==a2) continue;
					else {a4=ax; break;}
				}
				if(a1 && a4) {dihats.push_back(a1); dihats.push_back(a2); dihats.push_back(a3); dihats.push_back(a4); dihseq.push_back(dihats); cout <<a1->toString()<<"-"<<a2->toString()<<"-"<<a3->toString()<<"-"<<a4->toString()<<"\n";}
			}
		}
		rotatedPoses=std::vector<Pose>(1,this->extractPose());
		double DEL=(2*PI/prec);
		Eigen::Vector3d xN(1,0,0),yN(0,1,0),zN(0,0,1);
		for(int i=1;i<prec;i++)
		{
			this->rotateMolecule(xN*DEL);
			for(int j=1;j<prec;j++)
			{
				this->rotateMolecule(yN*DEL);
				for(int k=1;k<prec;k++)
				{
					this->rotateMolecule(zN*DEL);
					if(dihman)
					{
						std::vector<Pose> dihposes=this->getAllDihedralPosesFrom(dihseq,ff);
						for(Pose p : dihposes) {rotatedPoses.push_back(p);}
					}
					else rotatedPoses.push_back(this->extractPose());
				}
				this->rotateMolecule(zN*DEL);
				rotatedPoses.push_back(this->extractPose());
			}
			this->rotateMolecule(yN*DEL);
			rotatedPoses.push_back(this->extractPose());
		}
	}
public:
	void setLabel(const std::string& lab) {fragmentlabel=lab;}
	std::string getLabel() const {return fragmentlabel;}

	void alignByAtom(Atom* at,const Eigen::Vector3d& v)
	{
		Eigen::Vector3d trans=at->getPosition()-v;
		for(Atom* a : atoms) a->setPosition(a->getPosition()+trans);
	}
	/**@brief Rotate a fragment about a fixed atom position. The vector determines the direction, and it's magnitude the amount of rotation*/
	inline void rotateFragmentAboutAtom(const Eigen::Vector3d& n,Atom* centre) {this->rotateMolecule(n,centre->getPosition());}
	void alignFragmentToPose(Atom* o,Atom* s,Atom* t,const Eigen::Vector3d& ro,const Eigen::Vector3d& rs,const Eigen::Vector3d& rt)
	{
		alignByAtom(o,ro);
		Eigen::Vector3d ni=s->bondVectorTo(o).cross(s->bondVectorTo(t)),nf=Eigen::Vector3d(ro-rs).cross(Eigen::Vector3d(rt-rs));
		rotateFragmentAboutAtom(quickgeom::getRotation(ni,nf),o);
	}
	inline std::pair<Atom*,std::vector<Atom*>> getHandleAtom(Atom* ref) {return getHandleAtom(ref->toString());}
	std::pair<Atom*,std::vector<Atom*>> getHandleAtom(std::string ref="")
	{
		int hai; std::vector<int> sv;
		if(ref=="")
		{
			int ri=throwarandompoint(0,hands.size());
			hai=hands[ri];
			sv=supports[ri];
		}
		else
		{
			std::vector<int> al;
			std::vector<std::vector<int>> spts;
			for(int i=0;i<hands.size();i++)
			{
				if(atoms[hands[i]]->toString()==ref) {al.push_back(hands[i]); spts.push_back(supports[i]);}
			}
			if(al.size())
			{
				int ri=throwarandompoint(0,al.size());
				hai=al[ri];
				sv=spts[ri];
			}
			else return make_pair((Atom*)nullptr,std::vector<Atom*>());
		}
		std::vector<Atom*> alist; for(int& ind : sv) alist.push_back(atoms[ind]);
		return make_pair(atoms[hai],alist);
	}
	std::pair<Atom*,std::vector<Atom*>> getHandleAt(int i)
	{
		std::vector<Atom*> hsup; for(int& r : supports[i]) hsup.push_back(atoms[r]);
		return make_pair(atoms[i],hsup);
	}
	const std::vector<Pose>& getAllPoses() const {return rotatedPoses;}
	std::vector<Pose> produceRandomPoses(int num,Atom* cha,Eigen::Vector3d fixed,int devid=0)
	{
		std::vector<Pose> poses;
		int ai=indexOf(cha);
		while(poses.size()<num)
		{
			int posind=(int)RNG::throwarandompoint(devid,0,rotatedPoses.size());
			Pose rP=rotatedPoses[posind];
			Eigen::Vector3d delV=fixed-rP[ai];
			Pose fP; for(int i=0;i<rP.size();i++) fP.push_back(rP[i]+delV);
			poses.push_back(fP);
		}
		return poses;
	}
};

#ifdef USE_TORCH
#include "support/torchmol.h"

torchmol::AtomFeaturizer* atom_featurizer;
#ifdef SERIAL
static torchmol::GrowingDeNovoMolecule* molenvmt=nullptr; //(Molecule* root,int maxs,const AtomFeaturizer& atom_feat,int foc=-1)
#else
static std::vector<torchmol::GrowingDeNovoMolecule*> molenvmt;
#endif

#include "support/reinforce.hpp"
#ifdef SERIAL
static dnvreinforce::Reinforce* mypolicy;
#else
static dnvreinforce::EnsembleReinforce mypolicy;
#endif
#endif


//Defining Molecule class member functions here
Molecule::Molecule(Atom* a) : Molecule() {
	atoms.push_back(a);
	freeatoms=atoms;
	bonds.push_back(std::vector<int>());
	if(a->isHydrogen()) HC++;
	if(VARIEDCOUNT && VCLIST->contains(a)) HAC++;
}
Molecule::Molecule(std::vector<Atom*> ats) : Molecule() {
	atoms=std::vector<Atom*>();
	for(Atom* a : ats) atoms.push_back(new Atom(*a));
	for(int i=0;i<ats.size();i++) bonds.push_back(std::vector<int>());
	for(Atom* a : atoms)
	{
		if(a->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(a)) HAC++;
	}
}
Molecule::Molecule(std::vector<Atom*> ats,std::vector<std::vector<int>> bs)
{
	atoms=std::vector<Atom*>();
	for(Atom* a : ats) atoms.push_back(new Atom(*a));
	bonds=bs;
	freeatoms=std::vector<Atom*>();
	int K=0;
	for(Atom* a : atoms)
	{
		if(a->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(a)) HAC++;
		a->setValency(a->getStandardValency()-bonds[K].size());
		if(a->canBond()) freeatoms.push_back(a);
		K++;
	}
}

Molecule::Molecule(int n) : Molecule() {
	atoms=AtomSet(n);
	bonds=std::vector<std::vector<int>>(n);
}

Molecule::Molecule(const Molecule* o) : Molecule()
{
	for(int i=0;i<o->getSize();i++) atoms.push_back(new Atom(*(o->atoms[i])));
	for(Atom* a : atoms) {if(a->seek_valency()) freeatoms.push_back(a);}
	bonds=o->bonds;
	type=o->type;
	myBE=o->myBE;
	myHBE=o->myHBE;
	HC=o->HC;
	if(VARIEDCOUNT) HAC=o->HAC;
}

Molecule::Molecule(const char* file,bool eval)
{
	ifstream f;
	f.open(file,ios::in);//opening file containing the data set
	if(!f) {cout << "Requested file: "<<file<<"\n"; throw UnableToOpenFileException();}
	loadMoleculeFromTXT(f,eval,false);
	f.close();//file closed
}
Molecule::Molecule(std::istream& f,bool eval,bool cpt,const ForceField* ff) : Molecule() {loadMoleculeFromTXT(f,eval,cpt,ff);}
void Molecule::loadMoleculeFromTXT(std::istream& f,bool eval,bool cpt,const ForceField* ff)
{
	//std::cout.precision(std::numeric_limits<double>::max_digits10);
	if(!f) {cout <<"In Molecule::loadMoleculeFromTXT:\t"<<"Bad inputstream provided\n"; throw GeneralIOException(); }
	int N;
	f>>N;//first number in the file contains the number of atoms
	HC=0; //Hydrogen count
	if(VARIEDCOUNT) HAC=0;
	//atoms=new Atom*[N];
	atoms=std::vector<Atom*>();
	bonds=std::vector<std::vector<int>>();
	freeatoms=std::vector<Atom*>();
	orders=std::vector<std::vector<unsigned short int>>();
	for(int i=0;i<N;++i)//copying atom co-ordinates and other information from file
	{
		int v;
		double x,y,z,wt,c,s,e,r;
		std::string id;
		f>>id>>x>>y>>z>>v>>wt>>c>>s>>e>>r;
		//Atom* tat=new Atom(id,x/10.0,y/10.0,z/10.0,v,wt,c,s,e,r);
		#ifdef PLACEHOLDERMOLECULES
			Atom* tat=new Atom(id,x/10.0,y/10.0,z/10.0,v,wt,c,s,e,r);
		#else
			Atom* tat=(ff)?ff->getAtom(id,x/10.0,y/10.0,z/10.0):new Atom(id,x/10.0,y/10.0,z/10.0,v,wt,c,s,e,r);
		#endif
		//tat=new Atom(id,x/10.0,y/10.0,z/10.0,v,wt,c,s,e,r);
		//if(ff->getCategory("ring").contains(tat)) {tat->mustCycl=true; tat->setCyclized();}
		atoms.push_back(tat);
		bonds.push_back(std::vector<int>());
		if(tat->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(tat)) HAC++;
		//delete[] id;
	}
	std::string bonddat;
	getline(f,bonddat);
	if(eval)
	{
		if(cpt)
		{
			for(int i=0;i<N;i++)
			{
				getline(f,bonddat);
				//std::vector<int> bdat;
				std::stringstream ss(bonddat);
				int bn;
				while(!ss.eof())
				{
					ss >> bn;
					cout << bn << " ";
					if(!contains(bonds[i],bn)) makeBond(atoms[i],atoms[bn]);
					if(ss.eof()) break;
				}
				//bonds.push_back(bdat);
				cout << " as bonddat\n";
			}
		}
		#ifndef PLACEHOLDERMOLECULES
			for(int i=0;i<N;i++)
				for(int j=i+1;j<N;j++)
					if(atoms[i]->distanceFrom(atoms[j])<=(atoms[i]->seek_radius()+atoms[j]->seek_radius())/2.0 && atoms[i]->canBond() && atoms[j]->canBond()) makeBond(atoms[i],atoms[j]);
		#endif
    for(Atom* a : atoms) {if(a->seek_valency()) freeatoms.push_back(a);}
	}
	else
	{
		std::string bin;
		if(cpt)  {for(int i=0;i<atoms.size();i++) getline(f,bin);}
	}
	if(cpt)
	{
		std::string eners; getline(f,eners);
		cout << eners << " was read for energy\n";
		std::stringstream ens(eners);
		ens >> myBE >> myHBE;
	}
	this->describeStructure();
}
Molecule::Molecule(const std::string& infile,const ForceField& ff,const std::string& format,bool eval,bool warnf,char chain)
{
	ifstream inf; inf.open(infile,ios::in);
	if(!inf) {cout << "Requested file: "<<infile<<" was not found/could not be opened for reading\n"; throw FileNotFoundException();}
	loadFromInputStream(inf,ff,format,eval,warnf,chain);
	inf.close();
}
Molecule::Molecule(std::istream& infile,const ForceField& ff,const std::string& format,bool eval,bool warnf,char chain) {loadFromInputStream(infile,ff,format,eval,warnf,chain);}
void Molecule::loadFromInputStream(std::istream& inf,const ForceField& ff,const std::string& format,bool eval,bool warnf,char chain)
{
  if(format!="pdb" && format!="gro" && format!="mol2") {cout << "Error! Unknown format: "<<format<<"\n"; throw UnknownFileFormatException();}
	if(!inf) {cout <<"In Molecule::loadFromInputStream:\t"<<"Bad inputstream provided\n"; throw GeneralIOException(); }
  type=1;
  std::string line;
  int resn,atno=0;
  long line_no=0;
  atoms=std::vector<Atom*>();
  std::string atnm,rnam;
  double xx,yy,zz,chg;
  int N=0;
	HC=0;
	if(VARIEDCOUNT)
	{
		HAC=0;
		if(!VCLIST) {cerr << "ERR: "<<"Chosen varied counting but VCLIST is not set!\n"; throw VariedCountListNotSetException();}
	}
	char mode=(char)0;
	#ifndef BIG_MOLECULES
  int atomnumbers[100000];
	#else
	cout << "Big molecule mode!\n";
	int atomnumbers[1000000];
	int offset=0;
	#endif
  std::vector<int> connected;
	orders=std::vector<std::vector<unsigned short int>>();
  while(true)
  {
    line_no++;
    getline(inf,line);
    if(inf.eof()) break;
    if(eval && stringfx::indexOf("CONECT",line)==0)
    {
			line.erase(line.find_last_not_of(" \n\r\t")+1);
			if(!QUIET) cout << line << "\n";
      stringstream ss(line);
      ss >> atnm;
      int ind1,ind2; ss >> ind1;
      std::vector<int> blist;
      while(true)
      {
          ss >> ind2;
          blist.push_back(ind2);
          if(ss.tellg()==-1) break;
      }
      for(int i : blist)
			{
				if(atoms[atomnumbers[ind1]]->canBond() && atoms[atomnumbers[i]]->canBond() && !contains(connected,atomnumbers[i]))
				{
					makeBond(atoms[atomnumbers[ind1]],atoms[atomnumbers[i]]);
					if(!QUIET) cout << atoms[atomnumbers[ind1]]->toString() << " at index "<<atomnumbers[ind1]<<" bonds to "<<atoms[atomnumbers[i]]->toString()<< " at index "<<atomnumbers[i]<<"\n";
				}
			}
      connected.push_back(atomnumbers[ind1]);
    }
    if((format=="pdb" && stringfx::indexOf("ATOM",line)!=0 && stringfx::indexOf("HETATM",line)!=0) || (format=="gro" && (line[0]!=' ' && line[0]>'9' && line[0]<'0')) || (format=="mol2" && line[0]=='#')) continue;
    if(format=="pdb")
    {
      resn=std::stoi(line.substr(23,4));
      atno=std::stoi(line.substr(7,4));
      atnm=stringfx::trim(line.substr(12,4));
      rnam=stringfx::trim(line.substr(17,4));
      xx=std::stod(line.substr(31,7))/10;
      yy=std::stod(line.substr(39,7))/10;
      zz=std::stod(line.substr(47,7))/10;
    }
    else if(format=="gro")
    {
      if(line.length()<40 || line_no==1) continue;
      resn=std::stoi(line.substr(0,5));
      atnm=stringfx::trim(line.substr(10,5));
			atno=std::stoi(line.substr(15,5));
			#ifdef BIG_MOLECULES
			if(atno==99999)
			{
				atno+=offset;
				offset+=100000;
			}
			#endif
			rnam=stringfx::trim(line.substr(5,5));
			xx=std::stod(line.substr(21,8));
			yy=std::stod(line.substr(30,8));
			zz=std::stod(line.substr(37,8));
			if(!QUIET)
			{
				cout << resn <<" ";
	      cout << atnm <<" ";
	      cout << rnam <<" ";
	      cout << xx <<" ";
	      cout << yy <<" ";
	      cout << zz <<"\n";
			}
    }
		else if(format=="mol2")
		{
			cout <<"("<<mode<<") "<< line << "\n";
			if(line[0]=='@')
			{
				if(line.find("MOLECULE")!=std::string::npos) mode='M';
				else if(line.find("ATOM")!=std::string::npos) mode='A';
				else if(line.find("BOND")!=std::string::npos) mode='B';
				else mode='S';
				continue;
			}
			if(mode=='S' || !mode || line.length()<12) continue;
			switch(mode)
			{
				case 'M': continue;
				case 'A':
					atno=std::stoi(line.substr(3,4));
					atnm=stringfx::trim(line.substr(8,4));
					xx=std::stod(line.substr(18,8));
					yy=std::stod(line.substr(28,8));
					zz=std::stod(line.substr(38,8));
					rnam=stringfx::trim(line.substr(58,4));
					chg=std::stod(line.substr(69));
					break;
				case 'B':
					int ign;
					stringstream ss(line);
		      ss >> ign;
					std::string btp;
		      int ind1,ind2,ord=1; ss >> ind1 >> ind2 >> btp;
					if(atoms[atomnumbers[ind1]]->canBond() && atoms[atomnumbers[ind2]]->canBond())
					{
						try {ord=std::stoi(btp);} catch(std::exception& e) {}
						makeBond(atoms[atomnumbers[ind1]],atoms[atomnumbers[ind2]]);
						cout << "Orders: "<<orders.size() << "\n";
						orders[atomnumbers[ind1]].push_back(ord);
						orders[atomnumbers[ind2]].push_back(ord);
						if(!QUIET) cout << atoms[atomnumbers[ind1]]->toString() << " at index "<<atomnumbers[ind1]<<" bonds to "<<atoms[atomnumbers[ind2]]->toString()<< "\n";
					}
					continue;
			}
		}
    try
    {
      Residue rsd=ff.getResidue(rnam);
			std::pair<std::string,double> act=rsd.getActualWithCharge(atnm);
      atnm=get<0>(act);
			chg=get<1>(act);
			if(chg) mode='A';
			else {if(format!="mol2") mode=(char)0;}
    }
    catch(ResidueNotFoundException& ex)
		{
			if(!warnf)
				cout << "Warning: Residue '"<<resn<<"("<<rnam<<")': match not found. Assuming the molecule is not a regular protein/DNA molecule. (You will not see further warnings)\n";
			warnf=true;
			if(format!="mol2") mode=(char)0;
		}
		#ifndef PLACEHOLDERMOLECULES
	    Atom* match=nullptr;
	    for(Atom* a : ff.getAtomTypes())
	    {
	      if(stringfx::indexOf(a->toString(),atnm)!=0) continue;
	      if(!match || a->toString().length()>match->toString().length()) match=a;
	    }
	    if(!match)
			{
				cout << "Atom type match not found for: "<<line<<"for atom type:"<<atnm<<"\n";
				continue;
			}
	    atoms.push_back(ff.getAtom(match->toString(),xx,yy,zz));
		#else
			atoms.push_back(new Atom(atnm,xx,yy,zz,4));
		#endif
		if(mode=='A')
		{
			atoms[atoms.size()-1]->setCharge(chg);
			orders.push_back(std::vector<unsigned short int>());
		}
		if(atoms[atoms.size()-1]->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(atoms[atoms.size()-1])) HAC++;
    atoms[atoms.size()-1]->setResidue(rnam);
    atoms[atoms.size()-1]->setResidueNumber(resn);
    bonds.push_back(std::vector<int>());
    atomnumbers[atno]=N++;
    //if(contains(hotspotreses,resn)) {hotspotatoms.push_back(atoms[atoms.size()-1]); HN++;}
    //cout << line <<"\t"<<line.length()<<"\n";
  }
  //N=atoms.size();
	if(!QUIET) cout << N<<" atoms read into molecule under format: "<<format<<"\n";
  if(eval)
	{
    //cout << "Note: CONECT Records (if any) are omitted\n";
		for(int i=0;i<N;i++)
    {
      if(contains(connected,i)) continue;
			for(int j=i+1;j<N;j++)
				if(atoms[i]->distanceFrom(atoms[j])<=(atoms[i]->seek_radius()+atoms[j]->seek_radius())/2.0 && atoms[i]->canBond() && atoms[j]->canBond()) makeBond(atoms[i],atoms[j]);
    }
    for(Atom* a : atoms) {if(a->seek_valency()) freeatoms.push_back(a);}
	}
}

void Molecule::startTrials(int devid)
{
	if(trialmode) throw AlreadyInTrialModeException();
	oldatoms=atoms; oldbonds=bonds; oldfreeatoms=freeatoms; oldHC=HC;
	if(VARIEDCOUNT) oldHAC=HAC;
	oldbe=myBE; oldhbe=myHBE; olderrbe=errBE; olderrhbe=errHBE; oldcycles=cycles; valencies=std::vector<int>(); for(Atom* a : oldatoms) valencies.push_back(a->seek_valency());
	trialmode=true;
	
	#ifdef USE_TORCH
		#ifdef SERIAL
		molenvmt->startTrials();
		#else
		molenvmt[devid]->startTrials();
		#endif
	#endif
}
inline void Molecule::commit() {trialmode=false;}

bool Molecule::completeAttempt(ostream& logfile,const System& s,Atom* src,int trials,int nupper,int level,int levellim,bool restr,bool superseed,bool reseed,bool deseed,Molecule*& common_newroot,int common_seedlength,int devid)
{
	if(getEffectiveSize()>nupper) return false;
    if(level>levellim)
	{
		cout << "Size limit pre-trial reached\n";
		return false;
	}
    //if(!trialmode) startTrials();
    cout <<"For(ATOM): "<< src->toString() << "\n";
    logfile <<"For(ATOM): "<< src->toString() << "\n";
    Atom* nat=nullptr;
		try {nat=s.ff->selectAtomByRule(src,getBondedAtoms(src),nullptr,devid);}
		catch(RuleConstraintsPreventBondsException& ex) {std::cout << "Rule constraints prevent bonding!\n"; nat=nullptr;}
    if(!nat) {std::cout << "Warning: no atom found\n"; if(logfile) logfile << "No atom to bond to: "<<src->toString()<<"\n"; return false;}
		if(nat->seek_valency()>=2) {if(RNG::throwarandompoint(devid,0,1)<CYCPROB || (src->mustCycl && !src->isCyclized()) || (s.ff->hasCategory("ring") && s.ff->getCategory("ring").contains(nat))) nat->mustCycl=true;}
		//ADD Cyclization to atom-types here
    cout << "Trial: "<<nat->toString()<<"\t";
    logfile << "Trial: "<<nat->toString()<<"\t";
	#ifdef USE_TORCH
	float* rewptr=nullptr;
	#ifdef SERIAL
	if(mypolicy) rewptr=mypolicy->getRewardPointer();
	#else
	if(mypolicy[devid]) rewptr=mypolicy[devid]->getRewardPointer();
	#endif
	#endif
    Atom* nxt=nullptr;
    if(!attempt(s,src,nat,nxt,trials,restr,0,devid))
	{
		#ifdef USE_TORCH
		if(rewptr) *rewptr=ERR_REWARD;
		#endif
		delete nat;
		return false;
	}
    bondAtom(nat,src);

	#ifdef USE_TORCH
	if(rewptr && !nat->isHydrogen()) *rewptr=STEP_REWARD;
	#ifdef SERIAL
	molenvmt->updateLastAddedAtom();
	#else
	molenvmt[devid]->updateLastAddedAtom();
	#endif
	#endif
	
    if(nxt)
    {
      cout << nat->toString()<< "-"<<nxt->toString()<<"\n";
      logfile<<"CYCLIZED: " << nat->toString()<< "-"<<nxt->toString()<<"\n";
      makeBond(nxt,nat);
      /*cout << "Writing: cyclized.pdb\n";
      dumpMol("cyclized.pdb");*/
      int cycl1=indexOf(nat),cycl2=indexOf(nxt),indx;
      std::vector<int> psv=std::vector<int>(1,cycl1),ats; psv.push_back(cycl2);
      for(Atom* a : getBondedAtoms(nxt))
      {
        if(a==nat) continue;
        indx=indexOf(a);
        cout << indx << " is being tried for cycle\n";
        try
        {
          //cout << cycl1 <<","<<indx<<","<<psv<<"\n";
          ats=quickLocate(cycl1,indx,psv);
          cout << ats << " in the beginning\n";
          cout <<"Cycle("<<cycl1<<","<<cycl2<<")\n";
          logfile <<"Cycle("<<cycl1<<","<<cycl2<<") completed\n";
          for(int k =0;k<ats.size();k++)
					{
						atoms[ats[k]]->setCyclized(true);
						cout << atoms[ats[k]]->toString() <<"\t"; cout << atoms[ats[k]]->isCyclized() << ","; atoms[ats[k]]->cycl=true; cout << atoms[ats[k]]->isCyclized() <<"\n";
					}
          cout << ats << "in the end\n";
					std::vector<Atom*> psh; for(const int& ti : ats) psh.push_back(atoms[ti]);
					cycles.push_back(psh);
        }
        catch(NoPathPossibleException& ex) {}
      }
    }
    cout << "\n";
		if(superseed && getSize()==common_seedlength+1) common_newroot=new Molecule(this);
		#ifdef REALTIME
			this->dumpMol("realtime.pdb",s.ff);
			usleep(REALTIME);
		#endif
    while(!s.ff->isSatisfied(nat,getBondedAtoms(nat)))
    {
      if(!nat->canBond()) return false;
      if(!completeAttempt(logfile,s,nat,trials,nupper,level+1,levellim,restr,superseed,reseed,deseed,common_newroot,common_seedlength,devid)) return false;
    }
		if(s.ff->hasCategory("completion") && s.ff->getCategory("completion").contains(nat))
		{
    	while(nat->canBond()) {if(!completeAttempt(logfile,s,nat,trials,nupper,level+1,levellim,restr,superseed,reseed,deseed,common_newroot,common_seedlength,devid)) return false;}
		}
    describeStructure(logfile);
    describeStructure();
    /*startTrials();
    commit();*/
    cout << "Allowing commit of: "<<src->toString()<<"\n";
    return true;
}

Molecule::~Molecule() {for(Atom* a : atoms) delete a;}

inline std::pair<double,double> Molecule::calculateNonBondingEnergy(Atom* dummy,bool write)
{
	double eOT=0;
	if(write) cout << "Size: "<<this->getSize()<<"\n";
	for(int i=0;i<this->getSize();i++)
	{
		if(write)
		{
			cout << (long)(atoms[i]);
			cout << "\t'"<<atoms[i]->toString()<<"'\n";
		}
		eOT+=chemtools::getNonbondingPotential(dummy,atoms[i]);
	}
	return make_pair(eOT,0);
}

inline double required_energy(const Molecule*);
#ifndef MANUAL_SCORING
//Simplicity of energy calculation
inline double required_energy(const Molecule* m) {return m->binding_energy();}
#endif
//Added classes
/** \class Protein Molecule.hpp "graph/Molecule.hpp"
		This is the Protein class. DeNovo was initially designed to work mainly with protein receptors, and hence the class-name "Protein".<br/>
		The Protein class actually represents any standard macromolecule receptor on which the generation of binding ligands is desired. It is derived from the Molecule class, signfying that all Molecule functions are valid when called on a Protein as well.<br/>
		The generation procedure expects at-least one Protein object to be present in the System at the start of generation, and that this part of the system will remain unchanged throughout the run.<br/>
		The Protein class, in addition to storing atoms (similar to the Molecule class), also marks some atoms as "hotspot" atoms. The hotspot is the desired target domain of the protein on which the generation must be carried out. Hotspots can active sites of the macromolecule which are targetted for competetive inhibition, or some other allosteric sites.<br/>
		From DeNovo σ onwards, Proteins can be loaded directly from PDB or GRO files given that the atom-type names are correct
*/
class Protein : public Molecule
{
protected:
	std::vector<Atom*> hotspotatoms;
	int HN;
	double hsbe=0;
protected:
	Protein() : Molecule() {}
public:
	//Protein() : Molecule() {}
	/**@brief Recast a loaded molecule as a protein*/
	Protein(const Molecule* m) : Molecule(m) {type=1; HN=0;} //Protein
	/**@brief Load protein from TXT format (hotspot is unspecified). See also: Protein(const std::string& protein,const std::string& hotspot)*/
	Protein(const std::string& str) : Molecule(str.c_str()) {type=1; HN=0;} //Protein
	/**@brief Load the protein from TXT formats.
		 @details This function is to load the Protein from the old TXT format. One file has the entire protein (including the hotspot), and the other has only the hotspot atoms.<br/> See algorithm PDF for details of TXT format.
	*/
	Protein(const std::string& protfile,const std::string& hostspotfile) : Protein(protfile) {loadHotspot(hostspotfile);}
	/**@brief Load from a standard format file (PDB or GRO)
		 @details This function is used by the latest DeNovo versions. The input files are in standard formats<br/>This requires a ForceField to be passed separately.<br/> See Also: IncompleteForceField
		 @param[in] hotspotreses: Hotspot residue numbers. All atoms which belong to these residues are marked as "hotspot atoms"
		 @param[in] format: The file format which is being loaded. ("pdb" and "gro" are the supported values)
		 See Also: Molecule::Molecule(const std::string&,const ForceField& ff,const std::string&,bool,bool,char)
	*/
  Protein(const std::string& infile,const ForceField& ff,std::vector<int> hotspotreses,const std::string& format="pdb") : Molecule(infile,ff,format,false)
  {
		type=1;
		hotspotatoms=std::vector<Atom*>();
		for(Atom* a : atoms) {if(contains(hotspotreses,a->getResidueNumber())) hotspotatoms.push_back(a);}
		HN=hotspotatoms.size();
  }

	/**@brief Get the total number of atoms in hotspot residues*/
	virtual inline int getHotspotAtomCount() const {return HN;}

	/*void loadGrid(System& s)
	{
		//myGrid=tensor3D_int(wid,hei,dep); myGrid.fill_with_zero();
		cout << "Filling Protein grid\n";
		Eigen::Vector3i gp;
		Atom* ap;
		for(int i=0;i<getSize();i++) //Iterate through atoms of molecule
		{
			//cout << i << " "<< getSize()<<"\n";
			ap=this->getAtom(i);
			gp=castToGrid(this->getAtom(i)->getPosition(),s);
			//if(myGrid.get(gp(0),gp(1),gp(2))) cout << "Warning: Grid overlap at "<<gp(0)<<","<<gp(1)<<","<<gp(2)<<"\n";
			s.grid().modify(gp(0),gp(1),gp(2),1);
		}
		int tmp;
		cout << "Grid offset: "<<s.gridOffset(0)<<" "<<s.gridOffset(1)<<" "<<s.gridOffset(2) << "\n";
		cout << "Main protein done with "<<getSize()<<" atoms.\n";
		for(int i=0;i<HN;i++)
		{
			gp=castToGrid(this->hotspotatoms[i]->getPosition(),s);
			s.grid().modify(gp(0),gp(1),gp(2),2); //May need checking
		}
		cout << "Done. Protein and Hotspot are cast to grid.\n";
		cin >> tmp;
	}*/
	/**@brief Load the hotspot atoms from a TXT file
		 @details This is an old function (preserved for backward support). It is used to load the protein from DeNovo's TXT formar (see PDF)<br/>The hotspot atoms have to be separately placed in a TXT file which will be loaded.
	*/
	void loadHotspot(const std::string& hostspotfile)
	{
		ifstream f;
		f.open(hostspotfile,ios::in);//opening file containing the data set
		f>>HN;//first number in the file contains the number of atoms
    hotspotatoms=std::vector<Atom*>(HN);// hotspotatoms.reserve(HN);
		//hotspotatoms=new Atom*[HN];
		int i=0,v;
		double x,y,z,wt,c,s,e,r;
		for(;i<HN;++i)//copying atom co-ordinates and other information from file
		{
			char *id=new char[8];
			f>>id>>x>>y>>z>>v>>wt>>c>>s>>e>>r;
			hotspotatoms[i]=new Atom(id,x/10.0,y/10.0,z/10.0,v,wt,c,s,e,r);
		}
	}
	/* @brief Load the average deviatons from equilibrium for error calculation and more accurate calculation of potentials
		 @details The average deviation of each residue is loaded into each Atom (See Atom::posdev). This is used for calculating first-order deviation terms for error and adding second-order deviation terms to the value
	*/
	virtual void loadDevations(const std::vector<int>& resnums,const std::vector<double>& devs)
	{
		for(Atom* a : getAtoms())
		{
			for(int i=0;i<resnums.size();i++)
			{
				if(resnums[i]==a->getResidueNumber())
				{
					a->setPositionDeviation(devs[i]);
					break;
				}
			}
		}
	}

	/**@brief Compute the potential at a point in space*/
	virtual double getPotentialAt(const Eigen::Vector3d& pos) const
	{
		Atom* posat=new Atom("POS",pos(0),pos(1),pos(2),0,1,1,1,0,0.5);
		double pot=0;
		for(Atom* a : getAtoms())
			pot+=chemtools::getElectrostaticEnergy(a,posat);
		delete posat;
		return pot;
	}

	/**@brief Compute the potential at a point in space (3 individual parameters)*/
	inline virtual double getPotentialAt(const double& x,const double& y,const double& z) const {return getPotentialAt(Eigen::Vector3d(x,y,z));}

	/**@brief Get centre of geometry for hotspot*/
	Eigen::Vector3d getHotspotCentre() const
	{
		Eigen::Vector3d c(0,0,0);
		for(Atom* a : hotspotatoms)
			c+=a->getPosition();
		c/=hotspotatoms.size();
		return c;
	}

	/*Atom* generateSeed(const Eigen::Vector3d& O,const System& p,Atom* type=nullptr,double RC=1.0)
	{

	}*/
	/**@brief Generate seeds near the hotspot. DeNovo's default seeding algorithm
		 @details This seeding algorithm spreads seeds in the closest available free-space near the hotspot region of the protein. It keeps a minimum radius of spreading at 0.1nm and tries upto "RC". <br/> See the seeding algorithm in the PDF for more details
		 @param[in] p: The system. This system must contain this protein, and the temperature and ForceField of this system will be used. All atoms mentioned in seedtypes must be present in this ForceField
		 @param[in] seedtypes: The list of atom-types allowed to be used as seeds. For each position, an atom-type is selected randomly from the seed-list hence by repeating one atom-type, you can roughly double it's probability to appear as a seed (See seed reweighting in the PDF)<br/>As a general principle, it is useful to avoid hydrogen atoms as seeds (due to their small vanderwaal radius)
		 @param[in] maxseed: Maximum number of seeds to try to generate
		 @param[in] RC: The maximum radius upto which to try. If seeds are placed too far from the hotspot, it is likely that major contributions of interaction energy may not come from the hotspot region, hence the binders are likely to loose specificity (Default: 1nm)
		 @param[in] allowH: If set to true, hydrogen atoms are allowed to be seeds. It is set to true by default because it is assumed that the user will provide a valid seed-list and if hydrogen atoms are included in this list, they must be considered for seeds.
		 @param[in] dense: It determines the sparseness of distribution of the seeds. When dense mode is switched on, most of the seeds will be present in a very small (nearly spherical) region. Switching it off allows more spreading, but must be used with caution to ensure that the seeds do not leave the key regions (See PDF or source code for exact details - Default: off)
		 @param[in] densval: The rate of spreading (usually very small) when using "dense" mode. It allows finer control of the spreading during seeding.
	*/
	std::vector<Atom*> generateSeeds(const System& p,const std::vector<Atom*>& seedtypes=std::vector<Atom*>(),int maxseed=1,double RC=1.0,bool allowH=true,bool dense=false,double densval=0.005)
	{
		const ForceField& ff=*p.ff;
		Eigen::Vector3d cog=getHotspotCentre(); //Centre of geometry of the hotspot
		/*Eigen::Vector3d cog(0,0,0); //Centre of geometry
		for(int i=0;i<HN;i++)
			cog=cog+hotspotatoms[i]->getPosition();
		cog/=HN;*/
		Eigen::Vector3d fpos; double fen;
		fpos=cog;
		std::vector<Atom*> seeds;
		bool fnd=false;
		if(!QUIET) cout << "Using seed dispersion of "<<((dense)?densval:0.075)<<" per ring\n";
		for(double rc=0.1;rc<RC+0.1;rc+=(dense)?densval:0.075) //increments of 0.075A lengths (for radial check)
		{
			std::vector<Eigen::Vector3d> vecs;
			/*if(dense)*/ vecs=quickgeom::getRandomVectorsInShell(cog,rc,maxseed/2);
			//else vecs=quickgeom::getRandomVectorsInShell(cog,rc,maxseed/6);
			for(Eigen::Vector3d& v : vecs)
			{
				Atom* testa=(seedtypes.size())?new Atom(*randomSelect(seedtypes)):((!allowH)?ff.getRandomNonHAtomByType():ff.getRandomAtomByType());
				testa->setPosition(v);
				if(!p.strongWVCheck(testa)) {delete testa; continue;}
				if(get<0>(this->calculateNonBondingEnergy(testa))>0) {delete testa; continue;} // || boltzmannFactor(hsbe,p.temp)>throwarandompoint(0,1)) continue;
				if(p.ff->hasCategory("ring") && p.ff->getCategory("ring").contains(testa)) testa->mustCycl=true;
				seeds.push_back(new Atom(*testa));
				if(seeds.size()>=maxseed)
				{
					if(!QUIET)
					{
						for(Atom* a : seeds)
						cout << a->toString()<<"\t"<<a->getPosition().transpose()<<"\n";
					}
					return seeds;//{fnd=true; break;}
				}
				delete testa;
			}
			if(fnd) break;
		}
		if(!seeds.size())
		{
			cout << "Failed to generate seed in "<<RC<<" nm radius from: "<<cog(0)<<" "<<cog(1)<<" "<<cog(2)<<" - the COG (Centre of Geometry)\n";
			return std::vector<Atom*>();
		}
		/*for(auto& a : seeds)
			cout << get<0>(a)->toString()<<"\t"<<(get<1>(a)) << "\n";
		beta=1100;
		//for(int i=0;i<maxseed;i++) fseeds.push_back(get<0>(weighedSelectBoltzmann(seeds,p.temp)));
		beta=4.17;*/
		return seeds;
	}

	/**@brief Calculate the total non-bonding interaction energy of an atom with this protein
		 @details This method uses L-J and charge parameters of the supplied atom (and itself) to calculate the non-bonding interaction energy of this single atom with the Protein (or macromolecule).<br/> See DIEL (dielectric constant used in simulation). Refer to the "Charges" section in the algorithm document.
	*/
	virtual std::pair<double,double> calculateNonBondingEnergy(Atom* dummy,bool write=false) override
	{
		double r=get<0>(Molecule::calculateNonBondingEnergy(dummy,write));
		double thsbe=0;
		for(int i=0;i<HN;i++) thsbe+=chemtools::getNonbondingPotential(dummy,hotspotatoms[i]);
		return make_pair(r,thsbe);
	}
	/**@brief Calculate the total non-bonded interaction energy of a molecule to the Protein
		@details This method is the compiled energy DeNovo interaction calculation method. It calculates the interaction energy of the target Molecule object with itself (using L-J parameters and charges). It calculates the energy contribution of each atom and assigns it to the respective atoms. (See Atom::setEnergyContribution(double c))<br/> See also: DIEL  (dielectric constant used in simulation). For more details refer to the "Charges" section in the algorithm document.
	*/
	virtual std::pair<double,double> calculateNonBondingEnergy(Molecule* drug,bool write=false)
	{
		drug->myBE=0;
		drug->myHBE=0;
		for(Atom* a :drug->getAtoms())
		{
				const std::pair<double,double>& p=a->setEnergyContribution(calculateNonBondingEnergy(a,write));
		drug->myBE+=get<0>(p);
		drug->myHBE+=get<1>(p);
		}
		return make_pair(drug->myBE,drug->myHBE);
	}

	/**@brief Compute a grid with potential energy generated by the protein at each point in space
		@details Produces a grid of potential energy values for the protein. The grid is generated by sampling the protein at a regular grid of points.
	*/
	tensor3D_double* getPotentialGrid(int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
	{
		tensor3D_double* grid=new tensor3D_double(nx,ny,nz);
		double cx=xmin;
		double dx=(xmax-xmin)/nx,dy=(ymax-ymin)/ny,dz=(zmax-zmin)/nz;
		cout << "Bounds: "<<xmin<<" "<<xmax<<" "<<ymin<<" "<<ymax<<" "<<zmin<<" "<<zmax<<"\n";
		cout << "Steps:" << dx<<" "<<dy<<" "<<dz<<"\n";
		#ifndef SERIAL
		#pragma omp parallel for shared(grid)
		#endif
		for(int i=0;i<nx;i++)
		{
			double cy=ymin;
			for(int j=0;j<ny;j++)
			{
				double cz=zmin;
				for(int k=0;k<nz;k++)
				{
					std::cout << to_string(cx)+" "+to_string(cy)+" "+ to_string(cz) + " " + to_string(this->getPotentialAt(cx,cy,cz)) + "\n";
					grid->modify(i,j,k,this->getPotentialAt(cx,cy,cz));
					cz+=dz;
				}
				cy+=dy;
			}
			cx+=dx;
		}
		return grid;
	}
	/**@brief Compute the potential grid (default: Automatically compute the bounds)*/
	tensor3D_double* getPotentialGrid(int nx,int ny,int nz)
	{
		return getPotentialGrid(nx,ny,nz,this->getMin(0),this->getMax(0),this->getMin(1),this->getMax(1),this->getMin(2),this->getMax(2));
	}

	/**@brief \deprecated After calculating the non-bonded energy (see Protein::calculateNonBondingEnergy(Molecule*)), the hotspot energy can be got through this function*/
	inline double getHotspotBindingEnergy() {return hsbe;}
	/**@brief Get the entire list of atoms (as a list of pointers) which belong to the "hotspot" of this macromolecule (alternate method for a "const" Protein object)*/
  inline const std::vector<Atom*>& getHotspotAtoms() const {return hotspotatoms;}
	/**@brief Get the entire list of atoms (as a list of pointers) which belong to the "hotspot" of this macromolecule*/
  inline std::vector<Atom*> getHotspotAtoms() {return hotspotatoms;}

	/**@brief Generate Ligands (the main generation function to be called externally).
		 @details This function is used when you want to begin the generation from scratch (No seeds placed yet).<br/>
		 The remaining algorithm is described in generateLigands(const ForceField&,int,int,const std::string&,std::vector<Molecule*>,double,double,int,int,int,double,bool,bool,bool,bool,Molecule*,bool,int) and int the algorithm PDF
		 @param[in] seedtypes: The allowed atom-types to be used as seeds (std::vector<Atom*>). An atom-type is randomly picked from this list for each seed atom, hence if an atom type is present multiple times in the list, it has a higher chance of being selected (see seed-reweighting in algorithm PDF)
		 @param[in] seednum: The number of seeds to generate (int)
		 @param[in] seedfile: The file to which the initially generated seeds must be written (Default: seeds.pdb)
		 <br/>The other parameters are described in the other section (see above link)
	*/
	#ifndef FRAGGROW
	std::vector<System> generateLigands(const ForceField& ff,int n,int sizemin,const std::string& logfilename="dnv_protgen.log",std::vector<Atom*> seedtypes=std::vector<Atom*>(),int seednum=10,bool dense=false,double rmax=2,double ecut=1,int osccount=5,int tottries=250,int tries_per_lig=45,double temp=298,const std::string& seedfile="seeds.pdb",bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,const quickgeom::Container* contbox=nullptr,Molecule* initmol=nullptr,bool prune=false,float prunefrac=0.67)
	{
		if(!RNG::initialized) RNG::init();
		System tempS(ff); tempS.temp=temp;
		tempS.addMolecule(this);
		std::vector<Molecule*> rmols;
		std::vector<Atom*> seeds=generateSeeds(tempS,seedtypes,seednum,rmax,dense);
		std::vector<double> startbes,starthbes;
		for(int i=0;i<seeds.size();i++)
    {
			const auto& p=this->calculateNonBondingEnergy(seeds[i]);
      startbes.push_back(get<0>(p)); seeds[i]->setEnergyContribution(startbes[startbes.size()-1]);
      starthbes.push_back(get<1>(p)); seeds[i]->setHotspotEnergyContribution(starthbes[starthbes.size()-1]);
			Molecule* fmol=new Molecule(new Atom(*seeds[i]));
			fmol->myBE=startbes[startbes.size()-1];
			fmol->myHBE=startbes[starthbes.size()-1];
			rmols.push_back(fmol);
    }
		Molecule* seedm=new Molecule(0);
    for(Atom* a : seeds) seedm->addAtom(a);
    seedm->dumpMol(seedfile);
		return generateLigands(ff,n,sizemin,logfilename,rmols,rmax,ecut,osccount,tottries,tries_per_lig,temp,restr,superseed,reseed,deseed,contbox,initmol,prune,prunefrac);
	}
	#endif
	/**@brief The primary generation algorithm which is expected to be used externally. It has a large number of tunable parameters
		 @details The main generation begins from here. It uses seed molecules (which must be manually supplied), from which the growth continues.<br/>Each molecule accepted as monte-carlo step is written to indexed files under the name GENERATING_mol_*.pdb <br/>This method can also be used for <b>extensive pruning</b> (complete generation instead of quick filtering while pruning a molecule)
		 <br/>The generation is logged extensively and log files can go upto being as big as 2.5GB for a single long generation.
		 <br/>For the complete algorithm, see the PDF document
		 @param[in] ff: The forcefield to use for generation (See ForceField)
		 @param[in] n: The number of molecules to generate
		 @param[in] sizemin: The minimum size required. Growth will continue until this size is reached. If the valencies are satisfied eariler, the molecule is rejected unconditonally.
		 @param[in] logfilename: The name of the logfile to write to (Default: dnv_protgen.log). By default, in the parallel job, only the acceptances at each stage are recorded in the main log file. Other logfiles (one per each core used) are generated under the general name &lt;logfilename&gt;__&lt;core_number&gt;
		 @param[in] rootmols: The root molecules to begin with. These root molecules are divided among the cores (equally as far as possible), at the beginning of generation. If the number of cores assigned exceeds the seeds provided, each core is randomly assigned one eof the molecules ensuring that each molecule is assigned at-least once.
		 @param[in] rmax: Redundandant deprecated parameter. See generateLigands(const ForceField&,int,int,const std::string&,std::vector<Atom*>,int,bool,double rmax,double,int,int,int,double,const std::string&,bool,bool,bool,bool,Molecule*,bool,int)
		 @param[in] ecut: The energy tolerance (kJ/mol). If two molecules have energy within this tolerance level, the newer molecule is unconditonally accepted (irrespective of which is more stable). (Default 1kJ/mol)
		 @param[in] osccount: The number of oscillations after which a final batch of molecules is accepted. See algorithm PDF for more details.
		 @param[in] tottries: The total number of consecutive tries (for each core), before it stops the generation. Only if there are 'tottries' consecutive <b>failures in generation</b> (not rejections) does the generation stop. This is written in the log file of the respective core before quitting (Default: 250)
		 @param[in] tries_per_lig: The total number of tries per ligand in attempting to grow. (See algorithm flowchart in the PDF) (Default: 45)
		 @param[in] temp: Temperature at which the system is. Unit is Kelvin (Default: 298K)
		 @param[in] restr: Restrain the growth to stay within a limited distance from the "hotspot". This distance can be varied by changing the global variable RESTRAINDISTANCE. By default, RESTRAINDISTANCE=0.8 nm, but restrained growth is turned off.
		 @param[in] superseed: The superseed mode - Adds newer seeds during generation (based on interaction energies) to ensure faster convergence (See algorithm PDF)
		 @param[in] reseed: The reseed mode - Adds newer seeds during generation (based on interaction energies) to ensure faster convergence (See algorithm PDF)
		 @param[in] deseed: Deseed flag - When used with superseed or reseed, whenever a new seed is added, the source seed is removed. (See algorithm PDF)
		 @param[in] initmol: An inital template molecule to compete the generation with. If nothing is provided, the growth starts from scratch. Otherwise, the generation proceeds from this stage. This molecule is also taken as the template if the "prune" option is turned on
		 @param[in] prune: A flag that (if switched on) instructs DeNovo to prune the template molecule instead of growing DeNovo molecules to compete with it.
		 @param[in] prunedepth: The depth of pruning.
		 <br/>A sample invocation of this method is shown here (an excerpt from the standard generation CPP file)
		 @code{.cpp}
		 	ForceField ff("dnv/data/final_ff_parameters.ffin","dnv/data/categories.data");
			ff.loadBondParameters();
 			ff.loadRules("dnv/data/definitions.data");
		 	IncompleteForceField protff("dnv/data/protff.ffin"); protff.loadResidues("dnv/data/protein.rtp");
		 	std::vector<int> hsres={323,324,325,375,376,377,378,379,380,381,382,383,390,391,392,395,396,412,413,414,415,416,417,423,424,425,426,427,428,465,466,467,468,469,470,471,472,473,474,477,478,479,480,481,482};
			Protein p("protein_file.pdb",protff,hsres,format);
			std::vector<Atom*> seeds;
			seeds.push_back(ff.getAtom("CT1"));
			seeds.push_back(ff.getAtom("CT2"));
			seeds.push_back(ff.getAtom("CT3"));
			seeds.push_back(ff.getAtom("CA"));
			...

			...
			std::vector<System> mols=p.generateLigands(ff,num,size,"<log_name>.log",seeds,seednum,false,0.8,1,5,300,50,298,"seeds.pdb",<restrain?>,<superseed?>,<reseed?>,<deseed?>);
		 @endcode
		 The result is obtained as a collection of system into "mols", after which you can use System::getDrugMolecule() to get the generated drugs, and Molecule::dumpMol() to dump it to a PDB<br/>
		 Note: The relative paths are correct (relative to "dnv" wherever it was installed, but the path must be corrected for these lines for using standard DeNovo parameters in your CPP)
		 <br/>See Also: quickPrune(), Pruning - the procedure (in the algorithm PDF) for details on pruning
	*/
	#ifndef FRAGGROW
  std::vector<System> generateLigands(const ForceField& ff,int n,int sizemin,const std::string& logfilename="dnv_protgen.log",std::vector<Molecule*> rootmols=std::vector<Molecule*>(),double rmax=2,double ecut=1,int osccount=5,int tottries=250,int tries_per_lig=45,double temp=298,bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,const quickgeom::Container* contbox=nullptr,Molecule* initmol=nullptr,bool prune=false,float prunekeep=0.67)
	#else
	std::vector<System> generateLigands(const ForceField& ff,const FragmentSet& frags,int n,int sizemin,const std::string& logfilename="dnv_protgen.log",std::vector<Molecule*> rootmols=std::vector<Molecule*>(),double rmax=2,double ecut=1,int osccount=5,int tottries=MAX_FAILTRIES,int tries_per_lig=45,double temp=298,bool restr=false,bool superseed=false,bool reseed=false,bool deseed=false,const quickgeom::Container* contbox=nullptr,Molecule* initmol=nullptr,bool prune=false,float prunekeep=0.67)
	#endif
  {
		//Setup
		envprot=this;
		if(!RNG::initialized) RNG::init();
		#ifndef SERIAL
			int K[omp_get_max_threads()+1]; for(int i=0;i<=omp_get_max_threads();i++) K[i]=0;
		#else
			int K[2]; K[0]=0; K[1]=0;
		#endif
		if(!FFCHARGES)
		{
			netcomm::genff=&ff;
			netcomm::badsmileslog.open("badsmiles.log");
			netcomm::failures=&K[0];
			#ifndef NOZMQ
				#ifndef SERIAL
					zmqio::initThreads(cout,omp_get_max_threads());
				#else
					zmqio::initThreads(cout);
				#endif
			#endif
		}
		int tvar,maxsize=90,rootsize=-1;
    std::vector<System> ret;
		std::ofstream& commonlogfile=fileExists(logfilename);
		std::string cptfilename=logfilename;
		cptfilename=stringfx::replace(cptfilename,".log",".dnvcpt");
		if(std::ifstream(cptfilename))
		{
			cout << "Checkpoint "<<cptfilename<<" exists. Continue from here (Y/n)?: ";
			char myin; cin >> myin;
			if(myin!='N' && myin!='n') {ret=molfx::loadFromCheckpoint(cptfilename,ff); cout << ret.size() << " molecules loaded\n";}
			else cout << "Skipping (and overwriting) checkpoint file.\n";
		}

    System main(ff); main.temp=temp; main.addMolecule(this); if(contbox) main.setContainer(contbox);
		std::vector<System> coll;

		int osc=0;
		std::atomic<bool> lock(false);
		std::atomic<bool> lockcoll(false);

	commonlogfile << "Parameters:\n";
	commonlogfile << "Generating "<<n<<" drugs targetting a size of "<<sizemin<<((prune)?" above original size":"")<<" accumulating the best every "<<osccount<<" oscillations\n";
	commonlogfile << "Dielectric: "<<DIEL<<"\n";
	commonlogfile << "Restrained? ";
	if(restr) commonlogfile << "Y ("<<RESTRAINDISTANCE<<")\n";
	else commonlogfile << "N\n";
    commonlogfile << "Seeds\n";
    for(int i=0;i<rootmols.size();i++) {rootmols[i]->describeStructure(commonlogfile); commonlogfile <<"\twith energy: "<<rootmols[i]->myBE<<"("<<rootmols[i]->myHBE<<")\n";}
    commonlogfile.flush();
    //int rv;
		if(initmol)
		{
			System initsys=main;
			this->calculateNonBondingEnergy(initmol);
			initsys.addMolecule(initmol);
			initmol->type=0;
			initmol->dumpMol("inital_template.pdb");
			coll.push_back(initsys);
			if(prune) {rootmols=std::vector<Molecule*>(1,initmol); rootsize=initmol->getEffectiveSize(); maxsize=rootsize+sizemin+3;}
		}
		if(prune && !initmol && !rootmols.size()) prune=false; //Automatically disable pruning if an initial molecule is not given
		bool smilehandle=true;
		bool generate=true;
		if(!USE_PHYSICS) commonlogfile << "Primary Optimzation: On \n";
		else commonlogfile << "Primary Optimzation: Off \n";

		#ifndef SERIAL
		if(FORCE_CHARGES) commonlogfile << "Using uniform charge balancing to get integer net charges. \n";
		else commonlogfile << "WARN: Net charges per molecule may not be integers \n";

		if(!FFCHARGES) commonlogfile << "Default FF Charges: Off (using external charges) \n";
		else commonlogfile << "Default FF Charges: On \n";
		commonlogfile.flush();
		std::vector<Molecule*> threadseeds[omp_get_max_threads()];
		for(int i=0;i<omp_get_max_threads();i++) threadseeds[i]=std::vector<Molecule*>();
		for(int i=0;i<rootmols.size();i++) threadseeds[i%omp_get_max_threads()].push_back(rootmols[i]);
		for(int i=0;i<omp_get_max_threads();i++)  {if(!threadseeds[i].size()) threadseeds[i].push_back(randomSelect(rootmols));}
		commonlogfile << "Parallelized run begins with CPU count: "<<omp_get_num_procs()<<"\n";
		std::cout << "Parallelized run begins with CPU count: "<<omp_get_num_procs()<<", and ~"<<threadseeds[0].size()<<" seeds per core, with trial limit="<<tottries<<" per core\n";
		#ifndef NOINTERFACE
		std::cin >> tvar;
		#endif
		#pragma omp parallel
		{
			int id=omp_get_thread_num();
			cout << "Thread no "+to_string(id)+"\thas "+to_string(threadseeds[id].size())+" seeds\n";
			if(!FFCHARGES) cout << "Using external charges for thread "+to_string(id)+" "+to_string(generate)+" "+to_string(smilehandle)+"\n";
			else cout << "Using forcefield charges\n";
		#else
			int id=0;
			cout << "Serial Job\n";
		#endif

			int F=0,L=0,common_seedlength=0;
			bool dels,deltest=false;
			Molecule* testm=nullptr;
			System s1=System(ForceField());
			Molecule *chosenroot,*common_newroot=nullptr;
			double e1,e2,P,common_avgenergy;
			#ifndef SERIAL
			cout << "Starting on processor: "+to_string(id)+"\n";
			std::ofstream& logfile=fileExists(logfilename+"__"+to_string(id));
			logfile << "Starting on processor/thread: "+to_string(id)+" generation: "+to_string(generate)+", smiles: "+to_string(smilehandle)+"\n";
			logfile << "Using seeds: ";
			for(int i=0;i<threadseeds[id].size();i++) {threadseeds[id][i]->describeStructure(logfile); logfile <<"\twith energy: "<<threadseeds[id][i]->myBE<<"("<<threadseeds[id][i]->myHBE<<")\n";}
			logfile.flush();
			#else
				std::ofstream& logfile=commonlogfile;
			#endif

    	while(ret.size()<n)
	    {
				if(generate)
				{
					dels=false;
		      if(K[id]>tottries) {logfile << "All tries over. Quitting with results uptil now!\n"; commonlogfile <<id<<": All tries over. Quitting with results uptil now!\n"; logfile.close(); break;}
		      if(F>25) //Rejected 25 molecules in a row
		      {
		        //main.temp+=1;
		        commonlogfile << "25 Rejections (by energy) in a row. Heating system by 1K. Current temperature: "+to_string(main.temp)+"\n";
		        F=0;
		      }
		      s1=main;
					#ifndef SERIAL
						chosenroot=randomSelect(threadseeds[id]);
					#else
						chosenroot=randomSelect(rootmols);
					#endif
					deltest=false;
					testm=new Molecule(chosenroot);
					testm->type=0;
					#ifndef FRAGGROW
					if(prune)
					{
						std::cout << "PRUNE: "+std::to_string(prunekeep)+" as keep fraction\n";
						rootsize=testm->getEffectiveSize();
						testm->makeRandomPruneFragment(&ff,prunekeep,PRUNEKEEPHS,id);
					}
						//testm->makePruneFragment((int)RNG::throwarandompoint(id,0,prunedepth)+1);
					#endif
					#ifdef USE_TORCH
					#ifndef SERIAL
					molenvmt[id]->replaceMolecule(testm);
					#else
					molenvmt->replaceMolecule(testm);
					#endif
					#endif
					//logfile << "Loaded seed: "<<rootmols[rv]->toString()<<" with energy: "<<seeds[rv]->getEnergyContribution()<<"\n";

					s1.addMolecule(testm);
					if(superseed || reseed)
					{
						common_seedlength=testm->getSize();
						common_avgenergy=testm->myBE/common_seedlength;
					}
					if(!QUIET) cout << testm->getEffectiveSize() << " is effective size of testm to begin with\n";
					logfile << testm->getEffectiveSize() << " is effective size of testm to begin with\n";
							#ifndef FRAGGROW
		  			    if(prune) testm->samplegrow(logfile,s1,sizemin+rootsize,tries_per_lig,maxsize,LEVELLIMIT,restr,superseed,reseed,deseed,common_newroot,common_seedlength,id);
		  			    else testm->samplegrow(logfile,s1,sizemin,tries_per_lig,maxsize,LEVELLIMIT,restr,superseed,reseed,deseed,common_newroot,common_seedlength,id);
							#else
								testm->fragmentgrow(logfile,s1,frags,sizemin,tries_per_lig,maxsize,LEVELLIMIT,restr,id);
							#endif
					logfile << testm->getEffectiveSize() << " is effective size of testm after generation\n";
					bool accthis=true;
					if(!NOSEEDCHECK && !main.ff->isSatisfied(testm->getAtom(0),testm->getBondedAtoms(testm->getAtom(0)))) accthis=false;
					else if(!NOSEEDCHECK)
					{
						if(prune)
						{
							for(Atom* iat : testm->getAtoms())
							{
								if(!main.ff->isSatisfied(iat,testm->getBondedAtoms(iat)))
								{
									accthis=false;
									break;
								}
							}
						}
					}
					if(!accthis)
					{
						logfile << "Rejecting molecule due to unsatisfied seed atom\n";
						delete testm; deltest=true;
						continue;
					}
					if(!QUIET) cout << testm->getEffectiveSize() << " is effective size of testm after generation\n";

					if(superseed && common_newroot && testm->getSize()>common_seedlength)
					{
						if(!QUIET) cout << "Superseed: ";
						if(testm->getRandomBondableAtom() && testm->getAtom(common_seedlength)->getEnergyContribution()<common_avgenergy)
						{
							#ifndef SERIAL
								threadseeds[id].push_back(new Molecule(common_newroot));
								if(deseed) dels=true;
								if(!QUIET) cout << "Added new molecule of size = "<<common_newroot->getSize()<<"\n";
								common_newroot->dumpMol("newseed_"+to_string(id)+"_"+to_string(threadseeds[id].size())+".pdb");
							#else
								rootmols.push_back(new Molecule(common_newroot)); if(deseed) dels=true;
								if(!QUIET) cout << "Added new molecule of size = "<<common_newroot->getSize()<<"\n";
								common_newroot->dumpMol("newseed_serial_"+to_string(rootmols.size())+".pdb");
							#endif
						}
						else if(!QUIET) cout << "No change\n";
						delete common_newroot;
						common_newroot=nullptr;
					}
					if(reseed)
					{
						Atom *tatm,*fatm=nullptr;
						for(int i=common_seedlength;i<testm->getSize();i++)
						{
							tatm=testm->getAtom(i);
							if(tatm->isHydrogen()) continue;
							if(tatm->getEnergyContribution()<common_avgenergy) {if(!fatm || tatm->getEnergyContribution()<fatm->getEnergyContribution()) fatm=tatm;}
						}
						if(fatm) {rootmols.push_back(new Molecule(new Atom(*fatm))); commonlogfile << "Seed added: " << fatm->toString() << "\tTotal seeds: "<<rootmols.size()<<"\n";}
						if(deseed && fatm) dels=true;
					}
					if(dels) //Only activates with deseed
					{
						int delind=::indexOf(chosenroot,rootmols);
						delete rootmols[delind];
						//delete chosenroot;
						rootmols.erase(rootmols.begin()+delind);
						commonlogfile << "Seed erased.\n";
					}

					if(testm->getEffectiveSize()<=0) logfile << "It (testm) had negative size with "<<testm->HC<<" hydrogens\n"; //Check may be removed soon
					if(testm->getEffectiveSize()<sizemin)
					{
						cout << "Deleteing molecule\n"; logfile << "Deleteing molecule for short size: "<<testm->getEffectiveSize()<<"(H:"<<(testm->HC)<<") described by: ";
						testm->describeStructure(logfile);
						logfile<<", over required"<<sizemin<<"\n";
						delete testm; K[id]++; deltest=true;
						continue;
						//if(!smilehandle) continue;
					}
					K[id]=0;
					/*for(Atom* temp : testm->getAtoms())
					{
						switch(temp->toString()[0])
						{
							case 'H': case 'C': break;
							case 'N': testm->pE+=4.5; break;
							case 'O': testm->pE+=6; break;
							default: testm->pE+=10; break;
						}
						if(main.ff.getCategory("aromatic").contains(temp)) {testm->pE-=22;}
						else if(main.ff.getCategory("alkene").contains(temp)) {testm->pE-=5;}
						else if(temp->isCyclized()) {testm->pE-=12;}
					}*/ //Uncomment to penalize/reward by entropy (and polarity)
					if(!FFCHARGES)
					{
						ForceField nff(*main.ff);
						testm->calculateBondOrders3(nff);
						netcomm::processRequest(logfile,testm,id);
						//cout << "Completed generation of a molecule by: "+to_string(id)+"\n";
						logfile << "Completed generation of a molecule by: "+to_string(id)+"\n";
						//logfile << "Added request. Size now: "+to_string(netcomm::queueSize())+"\n";
						//while(testm!=netcomm::processRequest()) {}
					}
				}
				if(smilehandle)
				{
					//if(!testm) {lockcoll=false; continue;}
					if(!testm) continue;
					if(deltest) {cerr  << "ERR: Deleted molecule but still proceeded!\n"; cerr.flush();}
					if(USE_PHYSICS)
					{
						/*cout << "Completed rigid-body gradient descent\n";*/
						testm->dumpMol("GENERATING_mol_preopt."+to_string(id)+"_"+to_string(L)+".pdb",&ff);
						physics::gradientDescentPhysics(main,testm,RBDESCSTEPS);
						physics::atomisticGradientDescentPhysics(main,testm,FULLDESCSTEPS);
						logfile << "Minimized molecule by: "<<id<<"\n";
						logfile.flush();
					}
					if(FORCE_CHARGES) this->calculateNonBondingEnergy(testm);
					e1=required_energy(testm)+testm->pE;
					//if(e1!=e1 || e1==0.0000) {logfile << "Rejecting unstable molecule.\n"; delete testm; deltest=true; lockcoll=false; continue;}
					if(e1!=e1 || e1==0.0000) {logfile << "Rejecting unstable molecule.\n"; delete testm;  deltest=true; continue;}
					lock_guard<std::mutex> g1(myMutex1);
			    if(coll.size())
			    {
				  e2=required_energy(coll[coll.size()-1].getDrugMolecule())+coll[coll.size()-1].getDrugMolecule()->pE;
			      if(abs(e1-e2)<ecut) {coll.push_back(s1); logfile << "Similar energies ("<<e1<<","<<e2<<") encountered\n. "<<coll.size()<<" molecules collected.\n"; commonlogfile << to_string(id)+": Similar energies ("+to_string(e1)+","+to_string(e2)+") encountered\n. "+to_string(coll.size())+" molecules collected.\n"; lockcoll=false; continue;}
			      P=RNG::throwarandompoint(id,0,1);
				  if(P>boltzmannFactor(e1-e2,s1.temp)) {logfile << "Rejecting molecule of energy "<<e1<<" over more stable "<<e2<<" with P="<<boltzmannFactor(e1-e2,s1.temp)<<", throwing value in (0,1): "<<P<<"\n"; F++;  delete testm; deltest=true; lockcoll.store(false, std::memory_order_release); continue;}
			      else
			      {
			      	F=0;
			      	coll.push_back(s1);
					lockcoll.store(true, std::memory_order_release);
			      	//coll=std::vector<System>(1,s1);
			      	if(e1>e2) osc++;
			      	commonlogfile << "Accepting (possibly lower) energy molecule: "+to_string(e1)+" over "+to_string(e2)+"\n"; testm->describeStructure(commonlogfile);
			      	commonlogfile << "Oscillations upto now: "+to_string(osc)+"\nCommits upto now: "+to_string(ret.size())+"/"+to_string(n)+"\n";
					commonlogfile.flush();
			        testm->dumpMol("GENERATING_mol."+to_string(id)+"_"+to_string(++L)+".pdb",&ff);
			       }
			       //if(coll.size()>=n) {logfile << "Accumulated sufficient ligands: "<<n<<"\n"; return coll;}
			       if(osc>=osccount)
			       {
							lock_guard<std::mutex> g2(myMutex2);
							cout << "Oscillations reached cutoff="<<osccount<<". Will commit some molecules\n";
							commonlogfile << "Oscillations reached cutoff. Will commit some molecules\n";
							commonlogfile.flush();
							e1=0; double ex; int nc=0;
							for(System& s : coll)
							{
								 ex=required_energy(s.getDrugMolecule());
								 if(ex>0) continue;
								 e1+=ex;
								 nc++;
							}
							e1/=nc;
							commonlogfile << "Committing some molecules: Average binding energy (hotspot) = "+to_string(e1)+"\n";
							commonlogfile.flush();
			        std::vector<int> delinds;
			      	for(int i=1;i<coll.size();i++) //Purposefully skipping first entry as it will be a repeat from last generation
			        {
						//if(prune) {if(initmol->myBE<coll[i].getDrugMolecule()->myBE) continue;}
						if(required_energy(coll[i].getDrugMolecule())>e1 /*|| (initmol && initmol->myBE<coll[i].getDrugMolecule()->myBE)*/)
						{
							commonlogfile << "Deleteing collected molecule\n";
							coll[i].getDrugMolecule()->describeStructure();
							delinds.push_back(i);
							continue;
						}
						ret.push_back(coll[i]);
						//commonlogfile << "Committing with binding energy(array index: "<<i<<"): "<<required_energy(coll[i].getDrugMolecule())<<"\n";
						cout << "Committing with binding energy(array index: "<<i<<"): "<<required_energy(coll[i].getDrugMolecule())<<"\n";
			        }
			        for(int& i : delinds)
			        {
			          coll[i].getDrugMolecule()->describeStructure();
			        	delete coll[i].getDrugMolecule();
			        }
			        cout << "Current commit size: "<<ret.size()<<"\n";
			        if(ret.size()) coll=std::vector<System>(1,ret[ret.size()-1]);
			      	else coll=std::vector<System>();
							molfx::makeCheckpoint(cptfilename,ret,false);
			      	osc=0;
							lock.store(false, std::memory_order_release);
			    	}
			    }
			    else {coll.push_back(s1);}
				lockcoll.store(false, std::memory_order_release);
				}
	    }
		#ifndef SERIAL
		} //End the OMP parallel pragma
		#pragma omp barrier
		#ifndef NOZMQ
			zmqio::stopThreads(cout);
		#endif
		#endif
		if(!ret.size())
		{
			commonlogfile<<"ERR: Nothing was generated\n";
			cout << "WARN: No molecules could be generated.\n";
		}
		commonlogfile.close();
		return ret;
  }
	/**@brief Quickly prune a molecule in attempt to improve it's interaction to this protein by making modifications to the terminal atoms. See the pruning algorithm in the document for more details
		 @details Unlike a complete generation, this method quickly removes the terminal and regrows recursively for multiple times (as many iterations are specified), and selects the most stable alternative
		 @param[in] ff: The ForceField to use for this generation
		 @param[in] p: The molecule to prune
		 @param[in] logfilename: The name of the logfile into which the pruning is logged. This logging is not as extensive as generateLigands(). (Default filename: dnv_prune.log)
		 @param[in] iters: The number of iterations (and hence number of newly generated molecules) to try before concluding (Default: 100)
		 @param[in] expand: Allow expansion while pruning (i.e. allow the final molecule to be larger by 'expand' atoms than the initial molecule). Default: No expansion.
		 @param[in] tries_per_lig: See generateLigands()
		 @param[in] temp: The System temperature at which the growth is carried out (Default: 298K)
		 @param[in] prunedepth: Depth of pruning (See algorithm PDF or makePruneFragment())
	*/
	Molecule* quickPrune(const ForceField& ff,Molecule* p,std::string logfilename="dnv_prune.log",int iters=100,int expand=0,int tries_per_lig=45,double temp=298,int prunedepth=3,bool restr=false,bool forceexp=false)
	{
		envprot=this;
		int id=0;
		if(!p) return p;
		p->type=0;
		int tsize=p->getEffectiveSize()+expand;
		Molecule* curmol=p;
		this->calculateNonBondingEnergy(curmol);
		System main(ff); main.temp=temp; main.addMolecule(this);
		volatile bool lock=false;
		volatile int K=0;
		volatile int i=0;
		#ifndef SERIAL
			#pragma omp parallel private(id)
			{
				id=omp_get_thread_num();
				std::ofstream logfile; logfile.open(logfilename+"__"+to_string(id));
				//for(int i=0;i<iters/omp_get_max_threads()+1;i++)
		#else
			std::ofstream& logfile=fileExists(logfilename);
		#endif
		for(;i<iters;i++)
		{
			Molecule* tempm=new Molecule(curmol);
			tempm->makePruneFragment(prunedepth);
			//tempm->dumpMol("prunefrag.pdb");
			System s1=main;
			s1.addMolecule(tempm);
			tempm->samplegrow(logfile,s1,tsize,tries_per_lig,tsize+5,LEVELLIMIT,restr,false,false,false,NULLMOL,0,id);
			if(forceexp && tempm->getEffectiveSize()<tsize)
			{
				logfile << "Deleteing molecule due to small size "<<tempm->getEffectiveSize()<<" over "<<tsize<<"\n"; logfile.flush();
				delete tempm;
				continue;
			}
			if(!FFCHARGES)
			{
				ForceField nff=ff;
				tempm->calculateBondOrders3(nff);
				netcomm::processRequest(logfile,tempm,id);
		  }
			//if(forceexp && tempm->getEffectiveSize()<p->getEffectiveSize()+expand) {delete testm; continue;}
			if(USE_PHYSICS)
			{
				physics::gradientDescentPhysics(main,tempm,RBDESCSTEPS);
				physics::atomisticGradientDescentPhysics(main,tempm,FULLDESCSTEPS);
				this->calculateNonBondingEnergy(tempm);
			}
			if(tempm->myBE<curmol->myBE)
			{
				while(lock) {usleep(WAIT_TIME);}
				lock=true;
				if(tempm->myBE>curmol->myBE) {lock=false; delete tempm; continue;} //It is possible that while waiting, another thread commits a more stable molecule which may (unfortunately) be overwritten by this thread
				logfile << "Accepting molecule with energy: "+to_string(tempm->myBE)+" over less stable "+to_string(curmol->myBE)+"\n";
				curmol=tempm;
				tempm->dumpMol(logfile);
				lock=false;
				logfile.flush();
			}
			else
			{
				logfile << "Rejecting molecule with energy: "+to_string(tempm->myBE)+" over more stable "+to_string(curmol->myBE)+"\n";
				delete tempm;
			}
		}
		#ifndef SERIAL
			cout << "Process "+to_string(id)+" hit loop end\n";
			}
			#pragma omp barrier
		#endif
		return curmol;
	}
	#ifndef FRAGGROW
	std::vector<std::vector<System>> trajectoryWiseLigandGeneration(const ForceField& ff,const std::string& genname="dnv_trjgen",int trajcount=80,int base=10,int step=2,int maxsize=50,const std::vector<System>& chkpt=std::vector<System>(),std::vector<Atom*> defseeds=std::vector<Atom*>(),int defseednum=40,bool defdense=false,double defrmax=2,double defecut=1,int defosccount=5,int deftottries=250,int deftries_per_lig=45,double temp=298,int defpruneiters=150,int defprunedepth=-1,double prunescaleperstep=1/3.0)
	{
		if(defprunedepth==-1) defprunedepth=min(base/2,step*3);
		int prunedepth=defprunedepth;
		RNG::init();
		std::vector<std::vector<System>> trajes; //for(int i=0;i<trajcount;i++) trajes.push_back(std::vector<System>());
		std::ofstream logfile; logfile.open(genname+"_trj.log");
		std::vector<System> basemols;
		if(!chkpt.size())
		{
			basemols=generateLigands(ff,trajcount,base,genname+"_basegen.log",defseeds,defseednum,defdense,defrmax,defecut,defosccount,deftottries,deftries_per_lig,temp,genname+"_initseeds.pdb");
			molfx::makeCheckpoint(genname+"_base.cpt",basemols);
			logfile << "Baseline generation completed!\n";
		}
		else
		{
			basemols=chkpt;
			logfile << "Using provided checkpoint file as baseline generation with checkpoint base size:"+to_string(base)+"\n";
			int K=0;
			for(const System& s : chkpt)
			{
				s.getMolecules()[0]->type=0;
				this->calculateNonBondingEnergy(s.getDrugMolecule());
				//s.getDrugMolecule()->dumpMolAsTxt("startmol_"+to_string(K++)+".txt",true);
				s.getDrugMolecule()->dumpMol("startmol_"+to_string(K++)+".pdb",&ff);
			}
			//for(Atom* ap : basemols[0].getDrugMolecule()->getAtoms()) cout << ap->toString()<<"\t"<<ap->isHydrogen()<<" "<<basemols[0].getDrugMolecule()->indexOf(ap)<<"\t"<<basemols[0].getDrugMolecule()->getBonds()[basemols[0].getDrugMolecule()->indexOf(ap)]<<"\n";
		}
		logfile.flush();

		for(int i=0;i<basemols.size();i++) trajes.push_back(std::vector<System>(1,basemols[i]));
		int totsteps=(maxsize-base)/step;
		double prunedepthinc=0;
		System main(ff); main.addMolecule(this);
		for(int I=base+step;I<maxsize;I+=step)
		{
			logfile << "For step: "<<((I-base)/step)<<" of "<<totsteps<<" steps\n"; logfile.flush();
			for(int i=0;i<trajes.size();i++)
			{
				logfile <<"Beginning "<< (i+1) << " trajectory of "<<trajes.size()<<" total trajectories: "<<((double)(i*100)/trajes.size())<<"% complete\n";
				Molecule* curmol=trajes[i][trajes[i].size()-1].getDrugMolecule(); //Get last generated molecule
				if(!curmol) {logfile << "WARN: System "<<(i+1)<<" did not have (valid) base generation result at size: "<<I<<"\n"; logfile.flush(); continue;}
				curmol=this->quickPrune(ff,curmol,genname+"_prune_"+to_string(I)+".log__"+to_string(i),defpruneiters,step,deftries_per_lig,temp,prunedepth,false,true);
				System nsys=main.clone(); nsys.addMolecule(curmol);
				trajes[i].push_back(nsys);
			}
			logfile << "Completed step: "<<((I-base)/step)<<" of "<<totsteps<<" steps\n"; logfile.flush();
			std::vector<System> mid; for(int i=0;i<trajes.size();i++) mid.push_back(trajes[i][trajes[i].size()-1]);
			molfx::makeCheckpoint(genname+"_stage_"+to_string(I)+".cpt",mid);
			prunedepthinc+=prunescaleperstep*step;
			prunedepth=defprunedepth+prunedepthinc;
		}
		logfile.close();
		return trajes;
	}
	#endif

	/**@brief Get the radial charge distribution for this receptor
		 @details The charge distribution as a function of radius from a point can be calculated by this method.<br/>It will give the total charge for each distance. The resolution can be specified
		 @param[in] centre: The centre from which to calculate the Radial Charge Distribution
		 @param[in] resol: The resolution (in nm). The charge is summed in steps of "resol". (Default: 0.01nm)
		 @param[in] maxr: Maximum radius upto which to search. If the molecule is larger than this value, the program automatically expands to accomodate it. (Default 15nm)
		 @param[in] onlyhspt: Use only hotspot atoms for calculation of the charges (Default: off)
	*/
  std::vector<std::pair<double,double>> getRadialChargeDistribution(const Eigen::Vector3d& centre,double resol=0.01,double maxr=15,bool onlyhspt=false) //resol := resolution in nm
  {
    std::vector<std::pair<double,double>> ret;
    for(double i=resol;i<maxr;i+=resol) ret.push_back(make_pair(i,0));
    double d;
    int ind;
    std::vector<Atom*>& atlist=(onlyhspt)?hotspotatoms:atoms;
    for(Atom* a : atlist)
    {
      d=a->distanceFrom(centre);
      ind=(int)(d/resol);
      if(ind>ret.size())
			{
				cout << "Large protein: \n";
				return getRadialChargeDistribution(centre,resol,maxr+10);
			}
      get<1>(ret[ind])+=a->seek_charge();
    }
    return ret;
  }
};

//System class Methods
/*std::pair<double,double> System::calculateDeltaEnergy(Atom* dummy,Atom* src,const Bond& b,Molecule* srcM) const
{
	double dEB=srcM->calculateBondingEnergy(dummy,src,b);
	double dENB=calculateTotalNonBondingEnergy(dummy);
	cout << "Delta Energy Calculations complete\n";
	return make_pair(dEB,dENB);
}*/
/*double System::calculateTotalNonBondingEnergy(Molecule* srcM) const
{
	assert(srcM);
	double e=0;
	for(int i=0;i<srcM->getSize();i++) e+=calculateTotalNonBondingEnergy(srcM->getAtom(i),srcM);
	return e;
}
double System::calculateTotalNonBondingEnergy(Atom* dummy,Molecule* srcM=nullptr) const
{
	short int p;
	//cout << "Calculating NB energy of dummy atom: "<<dummy->toString()<<"\n";
	double dE=0,nbe=0,hsnbe=0;
	for(Molecule* m :mols)
	{
		nbe=m->calculateNonBondingEnergy(dummy);
		if(srcM)
		{
			if(m->type==1) //Protein
			{
				srcM->addBE(nbe); dummy->setEnergyContribution(nbe);
				hsnbe=((Protein*)m)->getHotspotBindingEnergy();
				srcM->addHBE(hsnbe); dummy->setHotspotEnergyContribution(nbe);
			}
		}
		dE+=nbe;
	}
	//cout << "Done\n";
	return dE;
}*/
void System::generateGrid(int marX=DEFAULT_MARX,int marY=5,int marZ=5,double res=10) //res is the resolution (10x resolution => 1nm has 10 lattice points)
{
	cout << "Begin grid construction\n";
	std::vector<double> xvals,yvals,zvals;
  for(Molecule* m : mols)
  {
  	for(int i=0;i<m->getSize();i++)
  	{
  		xvals.push_back(m->getAtom(i)->seek_x());
  		yvals.push_back(m->getAtom(i)->seek_y());
  		zvals.push_back(m->getAtom(i)->seek_z());
  	}
  }
	double mX,MX,mY,MY,mZ,MZ,wX,wY,wZ;
	mX=maths::min(xvals); MX=maths::max(xvals);
	mY=maths::min(yvals); MY=maths::max(yvals);
	mZ=maths::min(zvals); MZ=maths::max(zvals);
	cout <<"Will generate grid: xrange["<<mX<<":"<<MX<<"],yrange["<<mY<<":"<<MY<<"],zrange["<<mZ<<":"<<MZ<<"]\n";
	wX=(MX-mX)+2*marX; wY=(MY-mY)+2*marY; wZ=(MZ-mZ)+2*marZ;
	gridOffset=Eigen::Vector3d(mX-marX,mY-marY,mZ-marZ);
	int wid=(int)(wX*res+1),hei=(int)(wY*res+1),dep=(int)(wZ*res+1);
	myGrid=new tensor3D_int(wid,hei,dep);
	grid_width=wid; grid_height=hei; grid_depth=dep;
	grid_resolution=res;
	myGrid->fill_with_zero();
	cout << "Empty grid created.\n";
  for(Molecule* m : mols)
		m->loadGrid(*this);
	cout << "Grid filled with data\n";
}
inline System System::clone() const {return *this;}
System System::clone(std::vector<Molecule*> omit) const
{
	System ret=this->clone();
	ret.mols=std::vector<Molecule*>();
	for(int i=0;i<mols.size();i++)
	{
		if(contains(omit,mols[i]))
			ret.mols.push_back(new Molecule(mols[i]));
		else
			ret.mols.push_back(mols[i]);
	}
	return ret;
}
Eigen::Vector3i System::getClosestFreeSpace(int x,int y,int z,int maxtries=1e4)
{
	if(!occupied(x,y,z))
	{
		if(!occupied(x+1,y,z) && !occupied(x-1,y,z) && !occupied(x,y+1,z) && !occupied(x,y-1,z) && !occupied(x,y,z+1) && !occupied(x,y,z-1))
			return Eigen::Vector3i(x,y,z);
	}
	if(!maxtries) {cout <<"Sorry maxtries are up\n."; return Eigen::Vector3i(x,y,z);}
	int rd=throwarandompoint(0,6);
	switch(rd)
	{
		case 0:
			return getClosestFreeSpace(x-1,y,z,maxtries-1);
		case 1:
			return getClosestFreeSpace(x,y-1,z,maxtries-1);
		case 2:
			return getClosestFreeSpace(x,y,z-1,maxtries-1);
		case 3:
			return getClosestFreeSpace(x+1,y,z,maxtries-1);
		case 4:
			return getClosestFreeSpace(x,y+1,z,maxtries-1);
		case 5:
			return getClosestFreeSpace(x,y,z+1,maxtries-1);
		default:
			cout << "Warning: failed random int\n";
			return Eigen::Vector3i(x,y,z);
	}
}
void System::dumpSystem(const std::string& of) const
{
	ofstream f; f.open(of);
	for(Molecule* m : mols)
	{
		if(m->type==1) continue;
		m->dumpMol(f);
	}
	f.close();
	cout << "Dumped system to: "<<of<<"\n";
}
void System::dumpSystem(const std::string& of,const std::string& tf) const
{
	dumpSystem(of);
	ofstream f; f.open(tf,fstream::out);
	for(Molecule* m : mols)
	{
		if(m->type==1) continue;
		m->dumpMolAsTxt(f);
	}
	f.close();
	cout << "Dumped system(cpt) to: "<<tf<<"\n";
}
Molecule* System::getDrugMolecule() const
{
	for(Molecule* m : mols) {if(m->type!=1) return m;}
	return nullptr;
}
Molecule* System::popDrugMolecule()
{
	int ind=-1;
	for(int i=0;i<mols.size();i++)
	{
		if(mols[i]->type!=1)
		{
			ind=i;
			break;
		}
	}
	if(ind==-1) return nullptr;
	Molecule* ret=mols[ind];
	mols.erase(mols.begin()+ind);
	return ret;
}
Protein* System::getProteinMolecule() const
{
	for(Molecule* m : mols) {if(m->type==1) return (Protein*)m;}
	return nullptr;
}
bool System::strongWVCheck(Molecule* exm) const
{
	for(Molecule* m : mols)
	{
		if(m==exm) continue;
		for(int i=0;i<m->getSize();i++)
		{
			for(int j=0;j<exm->getSize();j++)
			{
				if(m->getAtom(i)->distanceFrom(exm->getAtom(j))<(m->getAtom(i)->seek_sigma()+exm->getAtom(j)->seek_sigma())/2) return false; //check fails: Atom clashes with protein
				//if(m->return_atom(i)->distanceFrom(exm->return_atom(j))<0.24) return false; //check fails: Atom clashes with protein
			}
		}
	}
	return true;
}
bool System::strongWVCheck(Atom* exa) const
{
	for(Molecule* m : mols)
	{
		//if(m==exm) continue;
		for(int i=0;i<m->getSize();i++)
		{
			if(m->getAtom(i)->isHydrogen()) continue;
			if(m->getAtom(i)->distanceFrom(exa)<(m->getAtom(i)->seek_sigma()+exa->seek_sigma())/2) return false; //check fails: Atom clashes with protein
			//if(m->return_atom(i)->distanceFrom(exa)<0.24) return false; //check fails: Atom clashes with protein
		}
	}
	return true;
}
void Molecule::dumpMol(ofstream& f,const ForceField* ff) const
{
	//if(type==1) return;
	Atom* a; double rnc;
	if(binding_energy()) f << "REMARK "<<binding_energy()<<"\t"<<hotspot_binding_energy()<<"\t"<<errBE<<"\t"<<errHBE<<"(error not impletemented yet -- Last two columns)\n";
	f << "REMARK VERSION "<<VERSION<<"\n";
	if(ff && (rnc=getNetCharge(*ff))) f << "REMARK Net charge: "<<rnc<<"\n";
	if(VARIEDCOUNT) f << "REMARK "<<HAC<<" 'heavy' atoms counted while running with varied counting scheme\n";
	#ifdef FRAGGROW
		for(MolecularFragment* mfg : componentFrags) f << "REMARK FRAGMENT FILE: "<<mfg->getLabel()<<"\n";
	#endif
	for(const auto& ring : cycles)
	{
		f << "REMARK CYCLE "<<ring.size()<<" with atoms: ";
		for(Atom* a : ring) f << indexOf(a) <<" ";
		f<< "\n";
	}
	const Category* halo=(ff && ff->hasCategory("halogen"))?&(ff->getCategory("halogen")):nullptr;
	//cout <<atoms.size()<<" atoms in molecule\n";
	//cout << (long)halo<<" as halogen category.\n";
	for(int i=0;i<atoms.size();i++)
	{
		a=atoms[i];
		if(!a) continue;
		//std::string tstr;
		char xc[10],yc[10],zc[10],tstr[100],ener[10],atind[20],resnum[10];
		std::string nform="%3.3f";
		sprintf(atind,"%5d",i+1); sprintf(resnum,"%4d",a->getResidueNumber()); sprintf(xc,nform.c_str(),a->seek_x()*10); sprintf(yc,nform.c_str(),a->seek_y()*10); sprintf(zc,nform.c_str(),a->seek_z()*10); if(atoms[i]->getEnergyContribution()<1e3) sprintf(ener,"%3.1f",atoms[i]->getEnergyContribution()); else sprintf(ener,"inf");
		//if(halo && halo->contains(a) && (a->toString()[0]=='C' || a->toString()[0]=='B')) sprintf(tstr,"ATOM%c  %4d %4s%4s A%4d    %8s%8s%8s       %5s          %c%c\n",((a->mustCycl)?' ':' '),i+1,a->toString().c_str(),a->getResidue().c_str(),a->getResidueNumber(),xc,yc,zc,ener,a->toString()[0],a->toString()[1]);
		//else sprintf(tstr,"ATOM%c %5d %4s%4s A%4d    %8s%8s%8s       %5s           %c\n",((a->mustCycl)?' ':' '),i+1,a->toString().c_str(),a->getResidue().c_str(),a->getResidueNumber(),xc,yc,zc,ener,a->toString()[0]);
		if(halo && halo->contains(a) && (a->toString()[0]=='C' || a->toString()[0]=='B')) sprintf(tstr,"ATOM%c %5s %4s%4s A%4s    %8s%8s%8s       %5s          %c%c\n",((a->mustCycl)?' ':' '),atind,a->toString().c_str(),a->getResidue().c_str(),resnum,xc,yc,zc,ener,a->toString()[0],a->toString()[1]);
		//sprintf(tstr,"ATOM%c %5s %4s%4s A%4s    %8s%8s%8s       %5s          %c%c\n",((a->mustCycl)?' ':' '),atind,a->toString().c_str(),a->getResidue().c_str(),a->getResidueNumber(),xc,yc,zc,ener,a->toString()[0],a->toString()[1]);
		else sprintf(tstr,"ATOM%c %5s %4s%4s A%4s    %8s%8s%8s       %5s          %c\n",((a->mustCycl)?' ':' '),atind,a->toString().c_str(),a->getResidue().c_str(),resnum,xc,yc,zc,ener,a->toString()[0]);
		if(!QUIET) cout << tstr;
		f << tstr;
		//f <<(a->toString())[0]<<" "<<a->seek_x()*10<<" "<<a->seek_y()*10<<" "<<a->seek_z()*10<<"\n";
	}
	for(int i=0;i<atoms.size();i++)
	{
		//auto v=atoms[i]->getBondedAtoms();
		if(!bonds[i].size()) continue;
		char anum[6];
		if(!QUIET) cout << "CONECT ";
		f << "CONECT ";
		sprintf(anum,"%4d",i+1);
		if(!QUIET) cout << anum;
		f << anum;
		for(int j=0;j<bonds[i].size();j++)
		{
			char conrec[20];
			sprintf(conrec," %4d",bonds[i][j]+1);
			if(!QUIET) cout << conrec;
			f << conrec;
		}
		if(!QUIET) cout << "\n";
		f<< "\n";
	}
}

void Molecule::rollback(const System& s,int devid)
{
  for(Atom* a : atoms)
	{
		if(!contains(oldatoms,a))
		{
			if(!QUIET) cout << "Deleting: "<<a->toString()<<"\n";
			delete a;
		}
	} //Freeing memory
  atoms=oldatoms;
  bonds=oldbonds;
  freeatoms=oldfreeatoms;
  for(int i=0;i<atoms.size();i++) atoms[i]->setValency(valencies[i]);
	myBE=oldbe; myHBE=oldhbe;
	errBE=olderrbe; errHBE=olderrhbe;
	HC=oldHC;
	if(VARIEDCOUNT) HAC=oldHAC;
	cycles=oldcycles;
	for(Atom* a : atoms)
	{
		if(!a->isCyclized() || a->isCycleStatic()) continue;
		bool fnd=false;
		for(const auto& cyc : cycles) {if(contains(cyc,a)) {fnd=true; break;}}
		if(!fnd) a->setCyclized(false);
	}
  /*double hse=0;
  for(Molecule* m : s.getMolecules())
  {
    if(m->type!=1) continue; //If not protein, ignore
    myBE=((Protein*)m)->calculateNonBondingEnergy(this);
    hse+=((Protein*)m)->getHotspotBindingEnergy();
  }
  myHBE=hse;*/
  #ifdef USE_TORCH
	#ifdef SERIAL
	molenvmt->rollback();
	#else
	molenvmt[devid]->rollback();
	#endif
  #endif
  trialmode=false;
}
inline double Molecule::getHotspotEnergy(const System& s,double v)
{
  if(v!=0) return v;
  for(Molecule* m : s.getMolecules())
  {
    v+=get<1>(((Protein*)m)->calculateNonBondingEnergy(this));
    //v+=((Protein*)m)->getHotspotBindingEnergy();
  }
  return v;
}
bool Molecule::attempt(const System& s,Atom* gP,Atom* nat,Atom*& nxt,int trials,bool restr,int nps,int devid)
{
	if(!QUIET) cout << trials << " ";
	bool acc=true;
	nxt=nullptr;
	Eigen::Vector3d pos;
	std::pair<Eigen::Vector3d,std::pair<double,double>> penv;
	try
	{
		if(s.ff->hasCategory("aromatic")) penv=chemtools::getAllowedPositions(gP,this,nat,s,20,10,restr,&(s.ff->getCategory("aromatic")),devid);
		else penv=chemtools::getAllowedPositions(gP,this,nat,s,20,10,restr,nullptr,devid);
		pos=get<0>(penv); nat->setPosition(pos);
	}
	catch(NoPositionAvailableException& ex) {acc=false; nps++;}
	if(!nat->isHydrogen() && acc) //If selected atom is accepted and not a hydrogen atom, trial to see if it is close enough for cyclization (Note that if it is, acceptance checks themselves ensure that it can cyclize)
	{
		for(Atom* ca : atoms)
		{
			if(ca==gP || ca==nat) continue;
			if(ca->isHydrogen()) continue;
			if(ca->distanceFrom(nat)<CYCCUT)
			{
				if(!nxt && ca->canBond() && nat->seek_valency()>1)  nxt=ca;
				else {acc=false; break;}
			}
		}
		//if(!acc) break;
	}
	if(acc)
	{
		if(get<0>(get<1>(penv)))
		{
			const auto& p=nat->setEnergyContribution(get<1>(penv));
			myBE+=get<0>(p); //nat->setEnergyContribution(get<0>(get<1>(penv)));
			myHBE+=get<1>(p); //nat->setHotspotEnergyContribution(get<1>(get<1>(penv)));
			/*if(nxt && s.getProteinMolecule()!=nullptr)
			{
				myBE+=nxt->setEnergyContribution(s.getProteinMolecule()->calculateNonBondingEnergy(nxt));
				myHBE+=nxt->setHotspotEnergyContribution(s.getProteinMolecule()->getHotspotBindingEnergy());
			}*/
			//myHBE+=getHotspotEnergy(s,get<1>(penv));
		}
		return true;
	}
	else
	{
		if(!QUIET && nps>=4) cout << "Trials found no acceptable position thrice. Exiting loop...\n";
		if(trials>0 && nps<4) {return attempt(s,gP,nat,nxt,trials-1,restr,nps+1,devid);} //Trial 'trial' number of times, but don't waste time if No position is available (thrice)
		else return false;
	}
}

static Molecule* merge(const Molecule* m1,const Molecule* m2)
{
  if(!m1) return new Molecule(m2);
  if(!m2) return new Molecule(m1);
  return new Molecule(m1);
}

//Completing chemtools
static std::vector<std::pair<Eigen::Vector3d,std::pair<double,std::pair<double,double>>>> chemtools::getAllAllowedPositions(Atom* src,Molecule* mol,Atom* natom,const System& sys,int num,int tries,bool restr,const Category* aromcat,int devid)
{
	//if(!mol) cout << "Null molecule pointer\n";
	//cout << "Failing at: "<<src->toString()<<"\n";
	mol->describeStructure();
	std::vector<Atom*> bondcol=mol->getBondedAtoms(src); //,angcol;
	double length=0,ang=0,tempen,proten=0;
  std::vector<double> dih;
  std::vector<Eigen::Vector3d> temp;
  bool angex=true;
  try {length=get<0>(sys.ff->getLengthParameters(src->toString(),natom->toString()));}
  catch(DataNotAvailableException& ex)
	{
		std::cout << "No data for bond "<<src->toString()<<"-"<<natom->toString()<<"\n";
		throw NoPositionAvailableException();
	}
	//if(!bondcol.size()) return make_pair(randomSelect(quickgeom::getRandomVectorsInShell(src->getPosition(),length,num*100)),make_pair(0.0,0.0));
  if(!bondcol.size()) {angex=false; temp=quickgeom::getRandomVectorsInShell(Eigen::Vector3d(0,0,0),length,num*20,devid);}
	std::vector<std::pair<Eigen::Vector3d,std::pair<double,std::pair<double,double>>>> ret;
  //std::vector<double> diheners;
	bool acc,ign=false,cyc,dihexist=false;
  double dihe,intlje=0,hsben=0;
  double vdis,ct,md=100;
  std::vector<std::pair<Atom*,std::vector<Atom*>>> atgeom;
  std::vector<Atom*> ignore;
  for(Atom* a : bondcol)
  {
    ignore.push_back(a);
    std::vector<Atom*> gats;
    for(Atom* a2 : mol->getBondedAtoms(a))
    {
      if(a2==src) continue;
      gats.push_back(a2); ignore.push_back(a2);
    }
    atgeom.push_back(make_pair(a,gats));
  }
  std::pair<Atom*,std::vector<Atom*>> rs;
  Eigen::Vector3d fixedbond,defbond;
  if(angex)
  {
    rs=randomSelect(atgeom);
  	fixedbond=src->bondVectorTo(get<0>(rs));
    ang=sys.ff->getAngleParamtersFailsafe(get<0>(rs)->toString(),src->toString(),natom->toString())[0];
    if(!(get<1>(rs).size())) {dihexist=false;}
    else
    {
      Atom* rs2=randomSelect(get<1>(rs));
      //dih.push_back(sys.ff->getDihedralParametersFailsafe(rs2->toString(),get<0>(rs)->toString(),src->toString(),natom->toString())[0]);

				//if(!natom->isHydrogen() && aromcat && aromcat->contains(src) && !src->isCyclized())
			if(!natom->isHydrogen() && aromcat && aromcat->contains(src) && !src->isCyclized() && (!ALLOW_BIPHENYL || !rs2->isCyclized()))
			{
				if(aromcat->contains(rs2)) {cout << "Forcing 0 dihedral for aromatic atom\n"; dih=std::vector<double>(1,0.0);}
				else {cout << "Forcing 180 (degree) dihedral for aromatic atom wrt to non-aromatic atom\n"; dih=std::vector<double>(1,PI);}
			}
			else
			{
      	auto pars=sys.ff->getDihedralParametersFailsafe(rs2->toString(),get<0>(rs)->toString(),src->toString(),natom->toString());
      	//dih=std::vector<double>(1,pars[0]*pars[2]);
      	//cout<<"\t"<<dih[0]<<"\t"<< rs2->toString()<<","<<get<0>(rs)->toString()<<","<<src->toString()<<","<<natom->toString()<<"\n";
      	dih=ForceField::expandMultiplicity(pars[0],pars[2]);
			}
      //cout << dih <<":\tchosen for "<<rs2->toString()<<","<<get<0>(rs)->toString()<<","<<src->toString()<<","<<natom->toString()<<"\n";
      defbond=get<0>(rs)->bondVectorTo(rs2);
      dihexist=true;
    }
  }
	//cout << "Found parameters: Length="<<length<<", Angle="<<(ang*57.28)<<", Dihedral="<<dih<<"\n";
  //cout << "Required angle: "<<ang<<"\n";
	while(ret.size()<num)
	{
    //ang=rule->predictAngle(randomSelect(bondcol),src,natom);
    if(angex)
    {
      if(!dihexist) {/*cout << "No dihedral constraints\n";*/ temp=quickgeom::getConePositions(fixedbond,length,ang,num*10,devid);} //Checked that angle is in radians.
      else
      {
        //cout << "\ttemp = "<<"quickgeom::getArcPositions("<<fixedbond.transpose()<<","<<defbond.transpose()<<","<<length<<","<<ang<<","<<dih<<","<<num<<");\n";
        temp=quickgeom::getArcPositions(fixedbond,defbond,length,ang,dih,num,devid);
      }
    }
		for(auto& v : temp)
		{
			acc=true;
      cyc=false;
      dihe=0; intlje=0; hsben=0;
			natom->setPosition(v+src->getPosition());
      //cout << "Position: "<<natom->getPosition().transpose()<<"\n";
      if(angex && abs(chemtools::getAngle(natom,src,get<0>(rs))-ang)>THETACUT) {acc=false; break;} //The positions are very closely clustered so don't even try the rest.
      intlje=0;
      hsben=0;
      if(!natom->isHydrogen())
      {
        for(Molecule* m : sys.getMolecules())
        {
          if(m->type==1) continue; //calculated below
          for(Atom* c : m->getAtoms())
          {
            if(c==src) continue;
            vdis=c->distanceFrom(natom);
            if(vdis<RCUT) {acc=false; break;}
            //if(vdis<md) md=vdis;
            if(c->isHydrogen()) {if(!ign) intlje+=getNonbondingPotential(c,natom); continue;}
            try {ct=sys.ff->getBondLength(c,natom,true);} catch(DataNotAvailableException& ex) {ct=0;}
						//if(vdis<md) {acc=false; break;}
          	ign=contains(ignore,c);
            if(vdis<CYCCUT) //How about replacing with energy criteria?
            {
              //cout << "Too close\n";
							if(cyc)
							{
								if(!QUIET) cout << "RINGREJ\t"<<natom->toString()<<" with "<<c->toString()<<": Clashing Rings\n";
								acc=false; break;
							}
              if(!ct)
							{
								if(!QUIET) cout << "RINGREJ\t"<<natom->toString()<<" with "<<c->toString()<<": No bond defined\n";
								acc=false; break;
							}
              if(!sys.ff->acceptNewBond(src,bondcol,natom) && !ign) {acc=false; break;}
              if((!c->canBond() && !ign) || natom->seek_valency()<2) {acc=false; break;}

              if(!ign && abs(vdis-ct)>CYCLENVAR)
							{
								if(!QUIET) cout << "RINGREJ\t"<<natom->toString()<<" with "<<c->toString()<<": Bad bond-length\n";
								acc=false; break;
							}
              auto bats=mol->getBondedAtoms(c);
              if(!ign && (!sys.ff->isSatisfied(c,natom,bats) || !sys.ff->isSatisfied(natom,c,std::vector<Atom*>(1,src))))
							{
								if(!QUIET) cout << "RINGREJ\t"<<natom->toString()<<" with "<<c->toString()<<": Unsatisfied atom type\n";
								acc=false; break;
							}
              for(Atom* n : bats)
							{
								if(sys.ff->calculateAngleEnergy(n,c,natom,true,1.5*THETACUT)!=0)
								{
									if(!QUIET) cout << "RINGREJ\t"<<natom->toString()<<" with "<<c->toString()<<": Bad angle\n";
									acc=false; break;
								}
							}
							if(!acc) break;
              cyc=ign=true;
							if(!QUIET) cout << "RINGACC\t"<<natom->toString()<<" with "<<c->toString()<<"\n";
            }
            if(!ign) intlje+=getNonbondingPotential(c,natom);
          }
        }
      }
			else
			{
				for(Molecule* m : sys.getMolecules())
        {
					if(m->type==1) continue;
          for(Atom* c : m->getAtoms())
          {
						if(c==src) continue;
						if(contains(ignore,c)) continue;
						vdis=c->distanceFrom(natom);
						if(vdis<RCUT*0.67) {acc=false; break;}
					}
				}
			}
      if(!acc) break;
			proten=0;
      for(Molecule* m : sys.getMolecules())
      {
        if(m->type==1)
        {
          for(Atom* pa : m->getAtoms()) //Can restrain atom selection here to ensure only nearby atoms are selected
          {
            if(!natom->isHydrogen() && !pa->isHydrogen() && pa->distanceFrom(natom)<PROTCUT) {acc=false; break;}
            tempen=chemtools::getNonbondingPotential(pa,natom);
            //if(throwarandompoint(0,1)>boltzmannFactor(tempen,sys.temp)) {acc=false; break;}
            intlje+=tempen;
						proten+=tempen;
          }
					if(proten>SINGLEATOMENERLIM) acc=false;
					if(!acc) {break;}
					if(!natom->isHydrogen())
					{
						if(sys.getContainer()) acc=sys.getContainer()->contains(natom->getPosition());
						if(acc && restr)
						{
							acc=false;
							for(Atom* hsa : ((Protein*)m)->getHotspotAtoms()) {if(hsa->distanceFrom(natom)<RESTRAINDISTANCE) {acc=true; break;}}
						}
						if(!acc) {break;}
					}
          for(Atom* pa : ((Protein*)m)->getHotspotAtoms()) hsben+=chemtools::getNonbondingPotential(pa,natom);
        }
      }
      if(cyc) {dihe=0; intlje=0;}
			if(intlje>SINGLEATOMINTENERLIM) acc=false;
			if(!acc) continue;
      //dihe+=intlje;
			for(auto& p : atgeom)
			{
        Atom* a=get<0>(p);
				//cout << getAngle(a,src,natom) << "\t" << (sys.ff.getAngleParamtersFailsafe(a->toString(),src->toString(),natom->toString()))[0] << "\n";
				//if(abs(getAngle(a,src,natom)-(sys.ff.getAngleParamtersFailsafe(a->toString(),src->toString(),natom->toString()))[0])>(10*PI/180.0)) {cout << "R\\a-1,"; acc=false; break;}
        if(!dihexist || natom->isHydrogen()) if(sys.ff->calculateAngleEnergy(a,src,natom,true,THETACUT)!=0) {acc=false; break;}
				//angcol=mol->getBondedAtoms(a);
				/*for(Atom* a2 : get<1>(p))
				{
					if(a2==src) continue;
					//if(mol->calculateDihedralEnergy(a2,a,src,natom)<DIHEDENERCUT) {acc=false; break;}
					//if(ff.calculateDihedralEnergy(a2,a,src,natom,true,DIHEDTHETACUT)) {acc=false; break;}
          dihe+=sys.ff.calculateDihedralEnergy(a2,a,src,natom,false);
				}*/
				if(!acc) break;
			}
      //cout << "\n";
			if(acc) {/*if(cyc) return make_pair(v,make_pair(0,0));*/ ret.push_back(make_pair(v+src->getPosition(),make_pair(intlje,make_pair(proten,hsben))));}
			if(ret.size()>=num || cyc) break;
		}
		tries--;
		if(tries<=0) break; //{mol->dumpMolAsTxt("faileddrug.txt"); break;}
	}
	return ret;
}
static std::pair<Eigen::Vector3d,std::pair<double,double>> chemtools::getAllowedPositions(Atom* src,Molecule* mol,Atom* natom,const System& sys,int num,int tries,bool restr,const Category* aromcat,int devid)
{
	auto ret=chemtools::getAllAllowedPositions(src,mol,natom,sys,num,tries,restr,aromcat,devid);
  if(ret.size())
	{
		auto& elr=weighedSelectBoltzmannFirst(ret,sys.temp,devid);
		if(!QUIET)
		{
			cout << "\nTrialled\n";
			cout << "Internal: "+to_string(get<0>(get<1>(elr)))<<"\n";
		}
		return make_pair(get<0>(elr),get<1>(get<1>(elr)));
	}
  //if(ret.size()) {cout << "\nTrialled\n"; return (randomSelect(ret));}
  else
  {
		if(!QUIET) cout << "#FAILED: No valid position available!!\n";
    /*Molecule* dmol=new Molecule(mol);
    for(Eigen::Vector3d& v : temp) {cout << v.norm() << "\n"; v+=src->getPosition(); dmol->addAtom(sys.ff->getAtom(natom->toString(),v(0),v(1),v(2)));}
    dmol->dumpMol("choke.pdb");
    delete dmol;*/
    //exit(1);
    throw NoPositionAvailableException();
  }
}
static std::vector<std::pair<Pose,double>> chemtools::getAllAllowedFragmentPositions(Atom* src,Molecule* mol,MolecularFragment* frag,std::pair<Atom*,std::vector<Atom*>> contri,const System& sys,int num,int tries,bool restr,const Category* aromcat,int devid) //num is num per each atom trialled. So num**3 trials per residue
{
	std::vector<std::pair<Pose,double>> ret;
	Molecule* core=mol;
	int srcInd=mol->indexOf(src);
	//std::pair<Atom*,std::vector<Atom*>> contri=frag->getHandleAt(handle);
	Atom* trialler1=new Atom(*get<0>(contri));
	auto atposes1=getAllAllowedPositions(src,mol,trialler1,sys,num,tries,restr,aromcat,devid);
	for(auto& p1s : atposes1)
	{
		Molecule* temp=new Molecule(core);
		System ms1(*sys.ff);
		temp->bondAtom(new Atom(*trialler1),temp->getAtom(srcInd));
		Atom* src2=temp->getLastAtom();
		Atom* trialler2=new Atom(*(get<1>(contri)[0]));
		auto atposes2=getAllAllowedPositions(src2,temp,trialler2,sys,num,tries,restr,aromcat,devid);
		delete trialler2;
		delete temp;
	}
	delete trialler1;
	return ret;
}
static std::pair<Pose,std::pair<double,double>> chemtools::getAllowedFragmentPositions(Atom* src,Molecule* mol,std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> frag,const System& s,int num,int tries,bool restr,const Category* aromcat,int devid) // Complete rewrite of above function
{
	MolecularFragment* cfrag=get<0>(frag);
	Atom* trialler=new Atom(*(get<0>(get<1>(frag))));
	auto atposes1=getAllAllowedPositions(src,mol,trialler,s,num,tries,restr,aromcat,devid);
	delete trialler;
	std::vector<std::pair<Pose,std::pair<double,double>>> posens;
	Protein* sysprot=nullptr;
	Atom* anchor=get<0>(get<1>(frag));
	for(Molecule* m : s.getMolecules()) {if(m->type==1) {sysprot=(Protein*)m; break;}}
	//cout <<"\t *** "<< atposes1.size() << " trial atom positions available ***\n";
	for(auto& ps : atposes1)
	{
		Eigen::Vector3d fixV=get<0>(ps);
		std::vector<Pose> fposes=cfrag->produceRandomPoses(num*cfrag->getEffectiveSize(),get<0>(get<1>(frag)),fixV,devid); //(int num,Atom* cha,Eigen::Vector3d fixed,int devid=0)
		for(const Pose& pose : fposes)
		{
			bool accpose=true;
			double posener=0;
			cfrag->assignPose(pose);
			if(get<1>(get<1>(frag)).size())
			{
				double ang1=chemtools::getAngle(src,anchor,get<1>(get<1>(frag))[0]); //,dih1=chemtools::getDihedral(src,anchor,get<1>(get<1>(frag))[0],get<1>(get<1>(frag))[1]);
				double expang=s.ff->getAngleParamtersFailsafe(src,anchor,get<1>(get<1>(frag))[0])[0];
				if(::abs(ang1)<0.9 || abs(ang1-expang)>THETACUT) continue; //No sharp angles
				if(get<1>(get<1>(frag)).size()>1)
				{
					double dih1=chemtools::getDihedral(src,anchor,get<1>(get<1>(frag))[0],get<1>(get<1>(frag))[1]);
					auto expdihps=s.ff->getDihedralParametersFailsafe(src,anchor,get<1>(get<1>(frag))[0],get<1>(get<1>(frag))[1]);
					std::vector<double> expdihs=ForceField::expandMultiplicity(expdihps[0],expdihps[2]);
					bool pass=false;
					for(double expdih : expdihs)
					{
						if(abs(expdih-dih1)<DIHVAR)
						{
							pass=true;
							break;
						}
					}
					if(!pass) continue;
				}
			}
			std::vector<Atom*> immn=cfrag->getBondedAtoms(anchor);
			for(Atom* fat : immn)
			{
				double reqang=s.ff->getAngleParamtersFailsafe(src,anchor,fat)[0];
				if(::abs(chemtools::getAngle(src,anchor,fat)-reqang)>THETACUT*2) //Weaker constraints for fragments
				{
					accpose=false;
					break;
				}
				/*std::vector<Atom*> secn=cfrag->getBondedAtoms(fat);
				for(Atom* fat2 : secn)
				{
					if(fat2==anchor) continue;
					auto reqdihps=s.ff->getDihedralParametersFailsafe(src,anchor,fat,fat2);
					double calcdih=chemtools::getDihedral(src,anchor,fat,fat2);
					std::vector<double> reqdihs=ForceField::expandMultiplicity(reqdihps[0],reqdihps[1]);
					bool pass=false;
					for(double expdihv : reqdihs)
					{
						if(abs(expdihv-calcdih)<DIHVAR)
						{
							pass=true;
							break;
						}
					}
					if(!pass) {accpose=false; break;}
				}*/
				if(!accpose) break;
			}
			if(!accpose) continue;
			for(Atom* fat : cfrag->getAtoms())
			{
				for(Atom* mat : mol->getAtoms())
				{
					if(mat->distanceFrom(fat)<RCUT/1.333)
					{
						accpose=false;
						break;
					}
				}
				if(!accpose) break;
			}
			if(!accpose) continue;
			if(sysprot)
			{
				std::pair<double,double> enV=sysprot->calculateNonBondingEnergy(cfrag);
				if(get<0>(enV)>4) {accpose=false; continue;}
				else posens.push_back(make_pair(pose,enV));
			}
			else posens.push_back(make_pair(pose,make_pair(0,0)));
		}
	}

	if(posens.size())
	{
		std::pair<Pose,std::pair<double,double>> finpos=weighedSelectBoltzmannFirst(posens,s.temp,devid);
		cfrag->assignPose(get<0>(finpos));
		if(sysprot) get<1>(finpos)=sysprot->calculateNonBondingEnergy(cfrag);
		return finpos;
		//return make_pair(get<0>(finpos),make_pair(get<1>(finpos),0));*/
	}
  else
  {
		if(!QUIET) cout << "#FAILED: No valid position available!!\n";
    throw NoPositionAvailableException();
  }
}
bool Molecule::reduce(const System& s,int trials,bool forced,int devid) //Perfect?
{
	cout << "Reduction begins: "<<((forced)?"forced":"unforced")<<"\n";
	Atom *ha=s.ff->getAtom("H"),*nxt,*tmp,*aa;
	double div,ebe,ehbe;
	int h=0;
	Molecule* p = s.getProteinMolecule();
	int sz=atoms.size(); if(sz<=2) return false;
	for(int I=0;I<sz;I++)
	{
		while(atoms[I] && atoms[I]->seek_valency())
		{
			//cout << atoms[I]->toString() << " to be reduced\t";
			double blen=0.09;
			Atom* ta=(forced)?s.ff->selectRandomBondableAtomTerminal(atoms[I],ha):s.ff->selectAtomByRule(atoms[I],this->getBondedAtoms(atoms[I]),ha,devid,nullptr);
			if(!ta || !ta->isHydrogen())
			{
				if(ta) {cout << "In reduce: No hydrogen atom to bond to: "<<atoms[I]->toString() << ". Will switch to other terminals\n"; blen=s.ff->getBondLength(atoms[I],ta,true);}
				else { cout << "In reduce: No terminal atom to bond to: "<<atoms[I]->toString() << "\n"; break;}
			}
			Atom *nha=nullptr,*eca1=nullptr,*eca2=nullptr,*dia=nullptr;
			int nhai=-1,eca1i=-1,eca2i=-1,diai=-1;
			for(int i=0;i<bonds[I].size();i++)
			{
				if(bonds[bonds[I][i]].size()<=1) continue;
				nha=atoms[bonds[I][i]]; nhai=bonds[I][i]; break;
			}
			if(!nha) //If there is no "heavy" atom connected (rare), choose a random hydrogen atom. If even that is not present, reject reduction altogether for this atom (Lone atom reduction is not allowed)
			{
				if(!bonds[I].size()) break;
				nha=atoms[bonds[I][0]];
				nhai=bonds[I][0];
			}
			//cout << "Using: "<<nha->toString()<<" as the first neighbour\n";
			Eigen::Vector3d base=atoms[I]->bondVectorTo(nha);
			for(int i=0;i<bonds[I].size();i++) //Find 2 more atoms which are bonded to reduction target
			{
				if(bonds[I][i]==nhai) continue;
				if(!eca1) {eca1i=bonds[I][i]; eca1=atoms[eca1i];}
				else {eca2i=bonds[I][i]; eca2=atoms[eca2i]; break;}
			}
			/*cout << "Looking for other neighbours: "<<(eca1!=nullptr)<<" and "<<(eca2!=nullptr);
			if(eca1) cout <<"\n"<< eca1->toString() << "\t";
			if(eca2) cout << eca2->toString();
			cout << "\n";*/
			Eigen::Vector3d ref; double dih=0;
			if(bonds[nhai].size()<=1) ref=Eigen::Vector3d(0,0,1);
			else
			{
				for(int i=0;i<bonds[nhai].size();i++)
				{
					if(bonds[nhai][i]==I) continue;
					diai=bonds[nhai][i];
					dia=atoms[diai];
					break;
				}
				ref=nha->bondVectorTo(dia);
				if(eca1 && dia) {dih=chemtools::getDihedral(eca1,atoms[I],nha,dia); } // cout << "Picked dih="<<dih<<" from calculation "<<eca1->getPosition().transpose()<<", "<<atoms[I]->getPosition().transpose()<<", "<<nha->getPosition().transpose()<<" and "<<dia->getPosition().transpose()<<"\n";}
			}
			int m=1;
			switch(atoms[I]->getStandardValency())
			{
				case 4:
					dih+=(2*PI)/3;
					m=3;
					break;
				case 3:
					dih+=PI;
					m=2;
					break;
				default:
					m=1;
					break;
			}
			//cout <<"Under dihedral data: "<<dih<<" with multiplicity "<<m<<"\n";
			//cout <<"Under getArcPositions: "<<base.transpose()<<"\t"<<ref.transpose()<<":\t"<<blen<<","<<((atoms[I]->getStandardValency()==4)?109.5:120)*(PI/180.0)<<","<<ForceField::expandMultiplicity(dih,m)<<"\n";
			std::vector<Eigen::Vector3d> tposes=quickgeom::getArcPositions(base,ref,blen,((atoms[I]->getStandardValency()==4)?109.5:120)*(PI/180.0),ForceField::expandMultiplicity(dih,m),32,devid);
			bool acc=true;
			Eigen::Vector3d fpos(0,0,0);
			double angs=0,nangs;
			for(Eigen::Vector3d& v : tposes)
			{
				ta->setPosition(atoms[I]->getPosition()+v);
				acc=true;
				nangs=0;
				for(int i=0;i<bonds[I].size();i++)
				{
					if(atoms[bonds[I][i]]->distanceFrom(ta)<RCUT/1.5) {acc=false; break;}
					nangs+=chemtools::getAngle(ta,atoms[I],atoms[bonds[I][i]]);
					//cout << chemtools::getAngle(ta,atoms[I],atoms[bonds[I][i]]) <<" for angle: "<<ta->toString()<<"*-"<<atoms[I]->toString()<<"-"<<atoms[bonds[I][i]]->toString()<<",\t";
				}
				//cout << " as nangs against "<<angs<<"\n";
				if(!acc || nangs<angs) continue;
				angs=nangs;
				fpos=ta->getPosition();
			}
			if(!angs) {cout << "Failed to reduce: "<<atoms[I]->toString()<<" with std. valency: "+to_string(atoms[I]->getStandardValency())+"\n"; break;}
			ta->setPosition(fpos);
			bondAtom(ta,atoms[I]);
			if(p)
			{
				const auto& ptr=((Protein*)p)->calculateNonBondingEnergy(ta);
				ebe=get<0>(ptr);
				ehbe=get<1>(ptr);
				myBE+=(ta->setEnergyContribution(ebe));
				myHBE+=(ta->setHotspotEnergyContribution(ehbe));
			}
			if(!QUIET) cout << "Added a hydrogen atom. Left: "<<atoms[I]->seek_valency()<<"\n";
			h++;
		}
	}
	if(!QUIET) cout << h << " hydrogen atoms added.\n";
	if(ha) delete ha;
	cout << "Reduction ends: "<<((forced)?"forced":"unforced")<<"\n";
	for(Atom* a : atoms)
	{
		if(a->getPosition()[0]!=a->getPosition()[0])
		{
			cout << "Bad Hydrogen Atom Placement\n";
			this->dumpMol("badred.pdb",s.ff);
			throw NoPositionAvailableException(); //Might remove this later on
		}
	}
	return true;
}
/*bool Molecule::reduce(const System& s,int trials,bool forced,int devid) //Improved
{
	cout << "Reduction begins: "<<((forced)?"forced":"unforced")<<"\n";
	Atom *ha=s.ff->getAtom("H"),*nxt,*tmp,*aa;
	int h=0;
	std::vector<Atom*> chn,chn2;
	double div,ebe,ehbe,angs=0;
	bool nx=true,fx=false,alx=false;
	int sz=atoms.size();
	Molecule* p = s.getProteinMolecule();
	for(int I=0;I<sz;I++)
	{
		fx=false;
		Atom* a=atoms[I];
		if(!a || !a->canBond()) continue;
		while(a->seek_valency())
		{
		chn=getBondedAtoms(a);
		if(!QUIET)
		{
			cout << a->toString() << ": checked.\n";
			cout << "Taken atom\n";
		}
		Eigen::Vector3d bnd(1,0,0),fixed(0,0,0);
		if(chn.size())
		{
			nxt=randomSelect(chn);
			bnd=nxt->bondVectorFrom(a);
			chn2=getBondedAtoms(nxt);
			if(chn2.size()>1) {tmp=a; while(tmp==a) tmp=randomSelect(chn2); fx=true; fixed=tmp->bondVectorFrom(nxt);}
		}
		if(!QUIET)
		{
			if(!s.ff->isSatisfied(a,chn))
			{
				cout << "Unsatisfied atom type: "<<a->toString()<<"\n";
			}
			cout << "Chosen: "<<a->toString()<<"\t";
		}
		Atom* ta=(forced)?s.ff->selectRandomBondableAtom(a,ha):s.ff->selectAtomByRule(a,getBondedAtoms(a),ha);
		if(!ta || !ta->isHydrogen())
		{
			if(!QUIET) cout << "!H\n";
			if(ta) delete ta; break;
		}
		int m=6;
		aa=nullptr;
		switch(a->getStandardValency())
		{
			case 4:
				for(Atom* taa : chn) {if(taa->isHydrogen()) {aa=taa; break;}}
				if(aa && tmp) {div=chemtools::getDihedral(aa,a,nxt,tmp)+(2*PI)/3; m=3;}
				else {div=PI/3; m=6;}
				break;
			case 3:
				div=0;
				m=2;
				break;
			default:
				div=0;
				m=1;
		}
		//div=(2*PI)/m;
		nx=true;
		alx=false;
		angs=0;
		double nangs=0;
		Eigen::Vector3d fpos(0,0,0);
			std::vector<Eigen::Vector3d> tposes=quickgeom::getArcPositions(bnd,(fx)?fixed:quickgeom::getRandomVectorInShell(Eigen::Vector3d(0,0,0),1,devid),0.09,((a->getStandardValency()==4)?109.5:120)*(PI/180.0),ForceField::expandMultiplicity(div,m),30,devid);
			for(Eigen::Vector3d& v : tposes)
			{
				nx=true;
				ta->setPosition(a->getPosition()+v);
				nangs=0;
				for(Atom* x : chn) {if(ta->distanceFrom(x)<RCUT/1.5) {nx=false; break;} nangs+=chemtools::getAngle(ta,a,x);}
				if(!nx || nangs<angs) continue;
				else {angs=nangs; fpos=ta->getPosition(); alx=true;}
			}
			if(!alx && forced && aa) {fpos=aa->getPosition()+quickgeom::getRandomVectorInShell(fpos,0.1,devid); alx=true;}
			if(alx)
			{
				ta->setPosition(fpos);
				bondAtom(ta,a);
				//if(ta->seek_x()!=ta->seek_x()) cout << "Failed to position hydrogen atom joining to: "<<a->toString()<<"("<<a->seek_valency()<<") using bond with "<<nxt->toString()<<"("<<nxt->seek_valency()<<") with dihedral from "<<tmp->toString()<<"\n";
				if(p)
				{
					const auto& ptr=((Protein*)p)->calculateNonBondingEnergy(ta);
					ebe=get<0>(ptr);
					ehbe=get<1>(ptr);
					myBE+=(ta->setEnergyContribution(ebe));
					myHBE+=(ta->setHotspotEnergyContribution(ehbe));
				}
				ta=(forced)?s.ff->selectRandomBondableAtom(a,ha):s.ff->selectAtomByRule(a,getBondedAtoms(a),ha);
				if(!QUIET) cout << "Added a hydrogen atom. Left: "<<a->seek_valency()<<"\n";
				h++; //break;
			}
		if(!alx)
		{
			if(!QUIET) cout << "Could not reduce with required dihedral:"<<div<<"\t: "<<a->toString()<<"\n";
			break;
		}
		if(ta) delete ta;
		tmp=nullptr;
		fixed=Eigen::Vector3d(0,0,0);
		}
	}
	cout << "Reduced. "<<h<<" hydrogens added!\n";
	delete ha;
	return true;
}*/

inline bool chemtools::Rule::satisfies(Atom* s,Atom* n,Molecule* src) const {return chemtools::Rule::satisfies(s,n,src->getBondedAtoms(s));}
inline bool chemtools::Rule::satisfies(Atom* s,Molecule* src,const ForceField* ff) const {return chemtools::Rule::satisfies(s,src->getBondedAtoms(s),ff);}
bool chemtools::Rule::satisfies(Atom* s,Atom* n,const std::vector<Atom*>& bnd) const
{
  if(group[0][0]=='@') return true;
	if(s->mustCycl && !s->isCyclized() && n->getStandardValency()<2)
	{
		if(!QUIET) cout <<"REJCYCAPP: Cycle formation impossible. Rejecting "<<n->toString()<<" as cyclizable atom type\n";
		return false;
	}
  //No ring check (New atom is being added - may form ring)
  if(!contains(group,n->toString())) return true;
  int c=0;
  for(Atom* a : bnd)
  {
    if(contains(group,a->toString())) c++;
    if(c>=max) return false;
  }
  return true;
}
bool chemtools::Rule::satisfies(Atom* s,const std::vector<Atom*>& bnd,const ForceField* ff) const
{
  /*cout << s->toString() << "\t";
  for(Atom* a : bnd) cout <<a->toString() << " ";
  cout << "\n";*/
  //if(ff) cout << s->toString()<<"\t"<<ff->getCategory("ring").contains(s) << "\n";
  /*try {if(ff && ff->getCategory("ring").contains(s) && !s->isCyclized()) return false;}
  catch(NoSuchCategoryException& ex) {cout << "WARN: Category 'ring' not available!!!\n";}*/
  if(group[0][0]=='@') return true;
  int c=0;
  for(Atom* a : bnd)
  {
    if(contains(group,a->toString())) c++;
    if(c>max) return false;
  }
  return (c>=min);
}
inline bool chemtools::Rule::satisfiesStrong(Atom* s,Atom* n,Molecule* src) const {return chemtools::Rule::satisfiesStrong(s,n,src->getBondedAtoms(s));}
bool chemtools::Rule::satisfiesStrong(Atom* s,Atom* n,const std::vector<Atom*>& bnd) const
{
  if(group[0][0]=='@') return true;
  int c=0;
  for(Atom* a : bnd)
  {
    if(contains(group,a->toString())) c++;
    if(c>=max) return false;
  }
  if(contains(group,n->toString())) return true;
  else return (c>=min);
}
//Completing forcefield
inline Atom* ForceField::selectAtomByRule(Atom* a,Molecule* m,Atom* exa,int devid,Atom* (*selector)(const ForceField&,const std::vector<std::string>&,Atom*,const std::vector<Atom*>&,int)) const {return selectAtomByRule(a,m->getBondedAtoms(a),exa,devid,selector);}
inline std::vector<std::string> ForceField::listAtomsByRule(Atom* a,Molecule* m,Atom* exa) const {return listAtomsByRule(a,m->getBondedAtoms(a),exa);}
inline std::vector<std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>>> ForceField::listFragmentsByRule(Atom* src,Molecule* m,const FragmentSet& frags,Atom* exa,int devid) const {return listFragmentsByRule(src,m->getBondedAtoms(src),frags,exa,devid);}
std::vector<std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>>> ForceField::listFragmentsByRule(Atom* src,const std::vector<Atom*>& bnd,const FragmentSet& frags,Atom* exa,int devid) const
{
	std::vector<std::string> alats=listAtomsByRule(src,bnd,exa);
	std::vector<std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>>> ret;
	for(MolecularFragment* frag : frags)
	{
		std::vector<std::pair<Atom*,std::vector<Atom*>>> possib;
		for(std::string ta : alats)
		{
			auto pres=frag->getHandleAtom(ta);
			if(get<0>(pres)) possib.push_back(pres);
		}
		//cout << src->toString()<<" with choices in: {"; //<< possib.size()<< " as size\n";
		//for(Atom* a : alats) cout << a->toString()<< " ";
		//cout <<"} has "<<possib.size()<<" options for fragment with label '"<<frag->getLabel()<<"'\n";
		if(!possib.size()) continue;
		ret.push_back(make_pair(frag,randomSelect(possib,devid)));
	}
	cout <<" \t As selection: "<<ret.size()<<" \n";
	if(!ret.size()) {cerr << "No options to grow using fragments! Source atom: "<<src->toString()<<"\n";}
	return ret;
}
inline std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> ForceField::selectFragmentByRule(Atom* s,Molecule* m,const FragmentSet& fs,Atom* temp) const {return selectFragmentByRule(s,m->getBondedAtoms(s),fs,temp);}

#ifdef FRAGGROW
std::vector<Atom*> Molecule::bondFragment(MolecularFragment* nfrag,int ind,Atom* src)
{
	//Atom* nacpy=natm;
	if(!nfrag->getAtom(ind)->canBond() || !src->canBond()) return std::vector<Atom*>();
	int cP=this->getSize(),K=0;
	std::vector<Atom*> ret;
	for(Atom* ea : nfrag->getAtoms())
	{
		if(ea->isHydrogen()) HC++;
		if(VARIEDCOUNT && VCLIST->contains(ea)) HAC++;
		this->addAtom(new Atom(*ea));
		for(int ae : nfrag->bonds[K]) bonds[K+cP].push_back(ae+cP);
		K++;
	}
	int sI=indexOf(src);
	bonds[sI].push_back(cP+ind); bonds[cP+ind].push_back(sI);
	src->addBond(); atoms[ind+cP]->addBond();
	for(int i=0;i<nfrag->getSize();i++)
	{
		ret.push_back(atoms[cP+i]);
		if(atoms[cP+i]->canBond()) freeatoms.push_back(atoms[cP+i]);
	}
	if(!src->canBond()) freeatoms.erase(std::remove(freeatoms.begin(), freeatoms.end(), src), freeatoms.end());
	return ret;
}

bool Molecule::attemptFragment(const System& s,Atom* gP,std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> nfrag,Atom*& pivat,int trials,bool restr,int nps,int devid)
{
	if(!QUIET) cout << trials << " ";
	bool acc=true;
	pivat=get<0>(get<1>(nfrag));
	Pose pos;
	std::pair<Pose,std::pair<double,double>> penv;
	try
	{
		if(s.ff->hasCategory("aromatic")) penv=chemtools::getAllowedFragmentPositions(gP,this,nfrag,s,20,10,restr,&(s.ff->getCategory("aromatic")),devid);
		else penv=chemtools::getAllowedFragmentPositions(gP,this,nfrag,s,20,10,restr,nullptr,devid);
		pos=get<0>(penv); //nat->setPosition(pos);
	}
	catch(NoPositionAvailableException& ex) {cout << "No position was found for: '"<<get<0>(nfrag)->getLabel()<<"'\n";  acc=false; nps++;}
	if(acc)
	{
		get<0>(nfrag)->assignPose(pos);
		if(get<0>(get<1>(penv)))
		{
			myBE+=get<0>(nfrag)->myBE;
			myHBE+=get<0>(nfrag)->myHBE; //nat->setHotspotEnergyContribution(get<1>(get<1>(penv)));
			//myHBE+=getHotspotEnergy(s,get<1>(penv));
		}
		return true;
	}
	else
	{
		if(!QUIET && nps>=8) cout << "Trials found no acceptable position 7 times. Exiting loop...\n";
		if(trials>0 && nps<8) {return attemptFragment(s,gP,nfrag,pivat,trials-1,restr,nps+1,devid);} //Trial 'trial' number of times, but don't waste time if No position is available (7 times)
		else return false;
	}
	return false;
}

MolecularFragment* Molecule::completeAttemptFragment(ostream& logfile,const System& s,Atom* src,const FragmentSet& fs,int trials,int nupper,int level,int levellim,bool restr,int devid)
{
	if(getEffectiveSize()>nupper) return nullptr;
	if(level>levellim)
	{
		cout << "Size limit pre-trial reached\n";
		return nullptr;
	}
	//if(!trialmode) startTrials();
	cout <<"For(FRAG): "<< src->toString() << "\n";
	logfile <<"For(FRAG): "<< src->toString() << "\n";
	std::pair<MolecularFragment*,std::pair<Atom*,std::vector<Atom*>>> nfrag=make_pair(nullptr,make_pair(nullptr,std::vector<Atom*>()));
	try {nfrag=s.ff->selectFragmentByRule(src,getBondedAtoms(src),fs);}
	catch(RuleConstraintsPreventBondsException& ex) {std::cout << "Rule constraints prevent bonding!\n"; get<0>(nfrag)=nullptr;}
	if(!get<0>(nfrag)) {std::cout << "Warning: no fragment found\n"; if(logfile) logfile << "No fragment to bond to: "<<src->toString()<<"\n"; return nullptr;}
	else {std::cout << "Chose fragment with label: '"<<get<0>(nfrag)->getLabel()<<"'\n";}
	//cout << "Trial: "<<nat->toString()<<"\t";
	logfile << "Trial: "; (get<0>(nfrag))->describeStructure(logfile); logfile<<"\t";
	Atom* pivat=nullptr;
	if(!attemptFragment(s,src,nfrag,pivat,trials,restr,0,devid)) return nullptr;
	std::vector<Atom*> added=bondFragment(get<0>(nfrag),get<0>(nfrag)->indexOf(pivat),src);
	assert(pivat);
	//std::vector<Atom*> added;
	cout << "\n";
	#ifdef REALTIME
		this->dumpMol("realtime.pdb",s.ff);
		usleep(REALTIME);
	#endif
	if(s.ff->hasCategory("completion"))
	{
		const Category* compcat=&(s.ff->getCategory("completion"));
		for(Atom* nat : added)
		{
			if(compcat->contains(nat)) while(nat->canBond()) {if(!completeAttemptFragment(logfile,s,nat,fs,trials,nupper,level+1,levellim,restr,devid)) return nullptr;}
		}
	}
	describeStructure(logfile);
	describeStructure();
	/*startTrials();
	commit();*/
	cout << "Allowing commit of: "<<src->toString()<<"\n";
	return get<0>(nfrag);
}
#endif

namespace atomfx
{
	inline Atom* selectRandomAtom(Molecule* m) {return randomSelect(m->getAtoms());}
	Atom* selectRandomNonHydrogenAtom(Molecule* m) {Atom* a; while((a=selectRandomAtom(m))->isHydrogen()) continue; return a;}
	bool ALWAYSTRUE(Atom* a) {return true;}
	bool ISHYDROGEN(Atom* a) {return a->isHydrogen();}
	std::vector<Atom*> selectAtoms(Molecule* m,bool (*f)(Atom*),unsigned int lim=0,bool sfl=false,bool invcond=false)
	{
		std::vector<Atom*> ret;
		std::vector<Atom*> atlist=m->getAtoms();
		if(sfl) {atlist=schuffle(atlist);}
		for(Atom* ta : atlist)
		{
			if((invcond && !f(ta)) || (!invcond && f(ta))) ret.push_back(ta);
			if(lim && ret.size()>=lim) break;
		}
		return ret;
	}

	inline std::vector<Atom*> selectRandomAtoms(Molecule* m,int n) {return selectAtoms(m,ALWAYSTRUE,n,true);}
	inline std::vector<Atom*> selectRandomNonHydrogenAtoms(Molecule* m,int n) {return selectAtoms(m,ISHYDROGEN,n,true,true);}
}
namespace molfx
{
	/**@brief Make a checkpoint file that saves a list of molecules*/
	static void makeCheckpoint(std::ostream& os,const std::vector<System>& syses,bool includeprot)
	{
		os << "#DNV generated checkpoint file\n";
		std::vector<Molecule*> unmols; //Unique molecules (by pointers)
		std::vector<std::vector<int>> indexing;
		for(const System& s : syses)
		{
			std::vector<int> myind;
			for(Molecule* m : s.getMolecules())
			{
				if(!includeprot && m->type==1) continue;
				int ind=indexOf(m,unmols);
				if(ind==-1) {unmols.push_back(m); myind.push_back(unmols.size()-1);}
				else myind.push_back(ind);
			}
			indexing.push_back(myind);
		}
		os << unmols.size() << "\n";
		for(Molecule* m : unmols) m->dumpMolAsTxt(os,true);
		os << syses.size() << "\n";
		for(std::vector<int>& inds : indexing) {for(int& el : inds) os<<" "<< el; os << "\n";}
		os << "#Checkpoint ends\n";
	}
	static void makeCheckpoint(std::ostream& os,const std::vector<Molecule*>& unmols)
	{
		os << "#DNV generated checkpoint file\n";
		os << unmols.size() << "\n";
		for(Molecule* m : unmols) m->dumpMolAsTxt(os,true);
		os << unmols.size() << "\n";
		for(int i=0;i<unmols.size();i++) os<<" "<< i<<"\n";
		os << "#Checkpoint ends\n";
	}
	static std::vector<System> loadFromCheckpoint(const std::string& cptfile,const ForceField& ff)
	{
		std::ifstream f; f.open(cptfile);
		std::vector<Molecule*> allmols;
		std::string line;
		int temp;
		while(!f.eof())
		{
			getline(f,line);
			cout << line << "\n";
			if(line[0]=='#') continue;
			else break;
		}
		std::stringstream ss(line);
		int nmol; ss >> nmol;
		cout << "Loading "<<nmol<<" molecules\n";
		while(allmols.size()<nmol) allmols.push_back(new Molecule(f,true,true,&ff));
		cout << "Molecules loaded\n";
		std::vector<System> syses;
		//while(!f.eof()){getline(f,line); cout << line << "\n";}
		getline(f,line);
		cout << line << "\n";
		ss=stringstream(line);
		int snum; ss >> snum;
		System s(ff);
		for(int i=0;i<snum;i++)
		{
			System t=s.clone();
			std::string ns="";
			getline(f,ns);
			cout <<i<<"\t'"<< ns << "'\n";
			ss=stringstream(ns);
			int ind;
			while(!ss.eof())
			{
				ss >> ind;
				cout << ind << " of "<<allmols.size()<<" molecules\n";
				t.addMolecule(allmols[ind]);
			}
			syses.push_back(t);
		}
		cout << syses.size() << " systems returned\n";
		return syses;
	}
	static std::vector<Molecule*> loadManyGROMolecules(std::istream& inf,const ForceField& ff,bool eval,bool warnf)
	{
		std::string fstr="",dump,temp;
		int s=0;
		std::vector<Molecule*> ret;
		while(!inf.eof())
		{
			getline(inf,dump);
			cout << dump << "\n";
			inf >> s; getline(inf,dump);
			//if(!s) break;
			for(int i=0;i<s;i++)
			{
				getline(inf,temp);
				fstr+=temp+"\n";
			}
			if(inf.eof()) break;
			getline(inf,dump); //Box vectors
			std::stringstream ss(fstr);
			ret.push_back(new Molecule(ss,ff,"gro",eval,warnf));
			fstr="";
		}
		return ret;
	}
	/**@brief Calculate the RMSD between two molecules
		 @details It uses the deviation in position of each atom from both the molecule to calculate the RMSD value.<br/> It <b>DOES NOT</b> reorder the atoms. The atoms are loaded in the same order as in the file from which the molecule was loaded (or if it was manually constructed, the order is the same as the order of adding atoms to the molecule).<br/>It maps the heavy atoms index-by-index for calculating the RMSD.<br/>
		 @param[in] mol1: The first molecule
		 @param[in] mol2: The second molecule
		 @param[in] inch: Include hydrogen atoms? If set to true, even hydrogen atoms will be considered for calculation (Default: false)
	*/
	double RMSD(const Molecule* mol1,const Molecule* mol2,bool inch=false)
	{
		double sqdistsum=0,sqdist;
	  int n=0,K=0;
	  if(mol1->getEffectiveSize()!=mol2->getEffectiveSize()) {cout << "Not same number of heavy atoms for RMSD! Quitting...\n"; exit(1);}
		if(inch) {if(mol1->getSize()!=mol2->getSize()) {cout << "Not same number of atoms (hydrogen count doesn't match) for RMSD! Ignoring hydrogen atoms\n"; inch=false;}}
	  for(int i=0;i<mol1->getSize();i++)
	  {
	    if(mol1->getAtom(i)->isHydrogen() && !inch) continue;
	    sqdist=pow(mol1->getAtom(i)->distanceFrom(mol2->getAtom(K)),2);
	    sqdistsum+=sqdist;
	    //cout << sqdistsum << "\n";
	    //cout << mol1->getAtom(i)->toString() <<" -> "<<mol2->getAtom(K)->toString()<<"\t"<<sqdist<<"\n";
	    K++;
	    //if(inch) while(K<mol2->getSize() && mol2->getAtom(K)->isHydrogen()) K++;
	    if(K>=mol2->getSize()) {cout << "Molecule ended!\n"; break;}
	    n++;
	  }
	  sqdistsum/=n;
		return sqrt(sqdistsum);
	}
	/*@brief Similar to RMSD(const Molecule*,const Molecule*,bool) but calculates the RMSD separately for each residue and returns an array of doubles with the RMSD of each residue by number
		@param[in] maxres: It specifies the largest possible residue number in the (macro)molecule. The array will have maxres elements with all indices which do not have a corresponding residue set to 0.
	*/
	double* RMSDByResidue(const Molecule* mol1,const Molecule* mol2,int maxres=10000,bool inch=false)
	{
		if(mol1->getEffectiveSize()!=mol2->getEffectiveSize()) {cout << "Not same number of heavy atoms for RMSD! Quitting...\n"; exit(1);}
		if(inch) {if(mol1->getSize()!=mol2->getSize()) {cout << "Not same number of atoms (hydrogen count doesn't match) for RMSD! Ignoring hydrogen atoms\n"; inch=false;}}
		double* ret= new double[maxres]; for(int i=0;i<maxres;i++) ret[i]=0;
		int atcnt[maxres]; for(int i=0;i<maxres;i++) atcnt[i]=0;
		for(int i=0;i<mol1->getSize();i++)
	  {
			if(mol1->getAtom(i)->isHydrogen() && !inch) continue;
			double sqdistsum=0,sqdist;
			ret[mol1->getAtom(i)->getResidueNumber()]+=pow(mol1->getAtom(i)->distanceFrom(mol2->getAtom(i)),2);
			atcnt[mol1->getAtom(i)->getResidueNumber()]++;
		}
		for(int i=0;i<maxres;i++)
		{
			if(ret[i]!=0) ret[i]=sqrt(ret[i]/atcnt[i]);
		}
		cout << maxres << " max residues used.\n";
		return ret;
	}
	static double closestApproach(const Molecule* m1, const Molecule* m2)
	{
		if(!m1->getSize() || m2->getSize()) throw EmptyMoleculeObjectException();
		double d=m1->getAtoms()[0]->distanceFrom(m2->getAtoms()[0]);
		if(!d) return 0;
		double r;
		for(Atom* a : m2->getAtoms())
		{
			for(Atom* b : m1->getAtoms())
			{
				r=a->distanceFrom(b);
				if(r<d) d=r;
			}
		}
		return d;
	}
	#ifndef NOJSON
	//void assignChargesFromJSON(Json::Value& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen=nullptr)
	/**@brief Use a JSON Object to assign charges to the atoms of the molecule. This method is disabled if compiled without JSON support*/
	void assignChargesFromJSON(Json::Value& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen)
	{
		cout << "Assign charges called (Json)\n";
		cout <<atorder.size() << " heavy atoms ordered\n";
		int myInd=0,hInd=atorder.size(); double chg=0; std::string temp,name;
		for(Atom* a : atorder)
		{
			if(halogen && halogen->contains(a)) name=(a->toString()[0]=='C')?"Cl":((a->toString()[0]=='B')?"Br":std::string(1,a->toString()[0]));
			else name=std::string(1,a->toString()[0]);
			temp=data[name+to_string(myInd++)].asString();
			if(temp=="" || temp=="null") throw AtomNotInJSONException();
			a->setCharge(std::stod(temp.substr(0,7)));
			std::vector<Atom*> attbonds=mptr->getBondedAtoms(a);
			for(Atom* a2 : attbonds)
			{
				if(!a2->isHydrogen()) continue;
				temp=data["H"+to_string(hInd++)].asString();
				if(temp=="null" || temp=="") throw AtomNotInJSONException();
				a2->setCharge(std::stod(temp.substr(0,7)));
				cout << "Set charge of "<<a2->toString()<< " to "<<a2->seek_charge() << "\n";
			}
			if(a->getStandardValency()>attbonds.size()) hInd+=(a->getStandardValency()-attbonds.size());
		}
	}
	#endif
	//void assignChargesFromData(const std::vector<std::pair<std::string,double>>& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen=nullptr)
	/**@brief Use the data provided in this vector to assign charges to atoms of the Molecule
		 @details The data can be input using either comma separated strings or as a JSON String. See txtio::loadChargesFromText(), txtio::loadChargesFromJSONString(), and jsonio::readJSON() <br/> The specifics of the I/O formats can be found in the PDF.
		 @param[in] data: The vector representing the data mapping each atom uniquely to a charge value. See PDF for format details.
		 @param[in] atorder: The vector with pointers to each Atom in the <b>same order</b> as the data. This consist of only the havy atoms and the total number of heavy atoms in both data and this vector are to be exactly same.
		 @param[in] mptr: The pointer to the Molecule to which these atoms (in atorder) belong.
		 @param[in] halogen: The halogen category. See Category, and the algorithm PDF section of categories for more information.
	*/
	void assignChargesFromData(const std::vector<std::pair<std::string,double>>& data,std::vector<Atom*> atorder,Molecule* mptr,const Category* halogen)
	{
		cout << "Assign charges called (data)\n";
		cout <<atorder.size() << " heavy atoms ordered\n";
		int myInd=0,hInd=atorder.size(); double chg=0; std::string temp,name;
		for(int i=0;i<atorder.size();i++)
		{
			/*if(halogen && halogen->contains(a)) name=(a->toString()[0]=='C')?"Cl":((a->toString()[0]=='B')?"Br":std::string(1,a->toString()[0]));
			else name=std::string(1,a->toString()[0]);*/
			//temp=data[name+to_string(myInd++)].asString();
			//if(temp=="" || temp=="null") throw AtomNotInJSONException();
			atorder[i]->setCharge(get<1>(data[i]));
			std::vector<Atom*> attbonds=mptr->getBondedAtoms(atorder[i]);
			for(Atom* a2 : attbonds)
			{
				if(!a2->isHydrogen()) continue;
				//temp=data["H"+to_string(hInd++)].asString();
				//if(temp=="null" || temp=="") throw AtomNotInJSONException();
				if(hInd>=data.size()) throw UnmatchingSMILESAndStructurePairException();
				a2->setCharge(get<1>(data[hInd++]));
				//cout << "Set charge of "<<a2->toString()<< " to "<<a2->seek_charge() << "\n";
			}
			if(atorder[i]->getStandardValency()>attbonds.size()) hInd+=(atorder[i]->getStandardValency()-attbonds.size());
		}
	}
}
namespace FFfx
{
	/**@brief Filter out atom-types of a particular element from a list of many*/
	std::vector<std::string> filterByElement(const ForceField& ff,const std::vector<std::string>& in,const std::string& el)
	{
		std::vector<std::string> ret;
		for(const std::string& s : in)
		{
			if(ff.getElementName(s)==el) ret.push_back(s);
		}
		return ret;
	}
	bool iterativeAtomTypeAssignment(const ForceField& ff,Molecule* m,const std::vector<Atom*>& atl,const Category* halo=nullptr,bool checkGeom=true,Atom* src=nullptr,std::vector<Atom*> visited=std::vector<Atom*>(),bool passgeom=true)
	{
		if(!src) src=m->getAtom(0);
		cout << "Started with atom type:"<<src->toString()<<" at index "<<m->indexOf(src)<<" using geometry type:"<<checkGeom<<"\n";
		//cout << "Iteration begins: "<<src->toString()<<"\n";
		//if(contains(visited,src)) return true;
		std::vector<Atom*> visitedProper=visited;
		std::vector<Atom*> neigh=m->getBondedAtoms(src);
		std::vector<bool> vis;
		std::vector<Atom*> shortatl;
		visited.push_back(src);
		for(Atom* na : neigh)
		{
			if(contains(visited,na))
			{
				vis.push_back(true);
				if(!shortatl.size())
				{
					try{
						auto ab=ff.getAllowedBonds(na);
						for(const auto& aba : ab) shortatl.push_back(aba.getSecond(na,ff));
						/*cout << na->toString()<<" allows bonding to: ";
						for(Atom* ta : shortatl) cout << ta->toString()<<","; cout << "\n";*/
					}
					catch(NoBondsAvailableException& excp) {continue;}
				}
			}
			else vis.push_back(false);
		}
		bool srcH=(halo)?halo->contains(src):false,satT=false,succ=false,stdata=false;
		if(!shortatl.size()) {shortatl=atl; stdata=true;}
		for(Atom* refa : shortatl)
		{
			//cout << src->toString()<<"\n";
			//cout << ((long)refa) << "\t" << ((long)src) << "\n";
			if(!refa) continue;
			if(refa->getStandardValency()!=src->getStandardValency()) continue;
			if(refa->toString()[0]!=src->toString()[0]) continue;
			else {if(halo && ((halo->contains(refa) && !srcH) || (!halo->contains(refa) && srcH))) continue;}
			src->assignAtom(refa,false);
			cout << "Assigned: "<<(m->indexOf(src)+1)<<": atom type: "<<src->toString()<<"="<<refa->toString()<<" with "<<refa->getStandardValency()<<" as required std. valency\n";
			satT=true;
			for(int t=0;t<neigh.size();t++)
			{
				if(vis[t])
				{
					try {ff.getLengthParameters(neigh[t]->toString(),src->toString(),true);}
					catch(DataNotAvailableException& ex) {satT=false; break;}
					if(checkGeom)
					{
						for(int u=t+1;u<neigh.size();u++)
						{
							if(!vis[u]) continue;
							if(ff.calculateAngleEnergy(neigh[t],src,neigh[u],true,THETACUT*2))
							{
								satT=false;
								if(!QUIET) cout << src->toString()<<":\t"<<neigh[t]->toString()<<"-"<<src->toString()<<"-"<<neigh[u]->toString()<<" with energy: "<<ff.calculateAngleEnergy(neigh[t],src,neigh[u],true,THETACUT*2)<<" => Geometry failure\n";
								break;
							}
						}
					}
				}
				if(!satT) break;
			}
			if(!satT) continue;
			for(int ni=0;ni<neigh.size();ni++)
			{
				Atom* na=neigh[ni];
				if(!vis[ni])
				{
					if(!iterativeAtomTypeAssignment(ff,m,atl,halo,(passgeom)?checkGeom:(!checkGeom),na,visited)) {satT=false; break;}
				}
			}
			if(!satT) continue;
			if(ff.isSatisfied(src,neigh)) {succ=true; break;}
			else {cout << "\tNOT satisfied: "<<src->toString()<<" at index "<<(m->indexOf(src)+1)<<"\n"; /*m->dumpMol("failed.pdb",&ff);*/}
			/*{
				bool geomsat=true;
				for(int ni1=0;ni1<neigh.size();ni1++)
				{
					for(int ni2=ni1+1;ni2<neigh.size();ni2++)
					{
						if(ff.calculateAngleEnergy(neigh[ni1],src,neigh[ni2],true,THETACUT*2))
						{
							geomsat=false;
							if(!QUIET) cout << src->toString()<<":\t"<<neigh[ni1]->toString()<<"-"<<src->toString()<<"-"<<neigh[ni2]->toString()<<" with energy: "<<ff.calculateAngleEnergy(neigh[ni1],src,neigh[ni2],true,THETACUT*2)<<" => Geometry failure\n";
							break;
						}
					}
					if(!geomsat) break;
				}
				if(geomsat) {succ=true; break;}
			}
			else continue;*/
		}
		if(!stdata) {for(Atom* ta : shortatl) delete ta;}
		if(checkGeom && !succ) {cout << "Relaxing geometric constraints because no atoms found:"<<src->toString()<<" at index "<<m->indexOf(src)<<"!\n"; return iterativeAtomTypeAssignment(ff,m,atl,halo,false,src,visitedProper,false);}
		return succ;
	}
	void assignAtomTypes(const ForceField& ff,Molecule* m,std::vector<Atom*> atypelist=std::vector<Atom*>(),bool checkGeom=true,const std::vector<std::pair<std::string,std::vector<std::string>>>& refinements=std::vector<std::pair<std::string,std::vector<std::string>>>())
	{
		const Category* halo=(ff.hasCategory("halogen"))?&(ff.getCategory("halogen")):nullptr;
		if(!halo) cerr << "assignAtomTypes (Molecule.hpp): WARN: No \"halogen\" category in supplied ForceField. Chlorine and Bromine will not be identified!!!\n";
		if(!atypelist.size()) atypelist=ff.getAtomTypes();
		std::vector<Atom*> initlist;
		for(Atom* ra : m->getAtoms()) ra->setStandardValency(m->getBondedAtoms(ra).size());
		if(iterativeAtomTypeAssignment(ff,m,atypelist,halo,checkGeom,m->getAtoms()[0],initlist)) cout <<"Assigned atom types successfully!\n";
		else cout << "Could not find suitable atom-type names\n";
		cout << "Refining atom types...\t";
		for(Atom* a : m->getAtoms())
		{
			std::string anam=a->toString();
			for(const auto& p : refinements)
			{
				if(anam==get<0>(p))
				{
					for(const auto& nn : get<1>(p))
					{
						cout << nn << "\n";
						Atom* na=ff.getAtom(nn);
						a->assignAtom(na,false);
						if(ff.isSatisfied(a,m->getBondedAtoms(a))) break;
						delete na;
					}
					break;
				}
			}
		}
		cout << "done\n";
	}
}
//Completing netcomm
/**
	This namespace deals with the data being passed between threads when using external charges.<br/>
	In this case, the generating thread (the thread that completes a generation successfully) adds it's newly generated molecule (along with the System object in which it was generated) to a queue. <br/>
	The queue is processed molecule by molecule by another asynchronous thread, which clears the queue parallely as it is filled.
*/
namespace netcomm
{
	/**struct MolRequest Molecule.hpp "graph/Molecule.hpp"
		 This struct holds the data necessary to process one request from a thread.
		 A queue of MolRequest objects is maintained in netcomm. This struct contains the molecule System and the ID of the thread which generated it.
	*/
	struct MolRequest
	{
		System mol; /**<The System containing the generated molecule*/
		int src; /**<The source thread ID*/
		MolRequest() : mol(ForceField()) {src=0;}
		/**@brief The simple MolRequest object constructor with all necessary data provided*/
		MolRequest(System m,int s) : mol(m) {src=s;}
	};
	static volatile Queue<MolRequest> molq(QUEUE_LIMIT); /**<The Queue which contains all the requests placed (and not processed) until now*/

	/**@brief Get the current queue size*/
	inline int queueSize() {return molq.getSize();}
	/**@brief Add a new request to the queue*/
	inline void addRequest(System m,int id) {while(queueSize()>=QUEUE_LIMIT-1) {cerr << "Error. Queue full, so waiting."+to_string(id)+"\n"; cerr.flush(); usleep(3000);} molq.push(MolRequest(m,id)); cerr << "Added request. Size now: "+to_string(queueSize())+"\n";}
	inline std::vector<MolRequest> drainQueue() {return molq.drain();}

	inline void processRequest(std::ostream& logfile,Molecule*& testm,int src)
	{
		//if(molq.isEmpty()) {logfile << "No data\n"; usleep(10000); return System(*genff);}
		//MolRequest req=molq.next();
		//Molecule* testm=req.mol.getDrugMolecule();
		if(!testm) {cerr << "WARN: No molecule in system\n"; return;}
		auto res=testm->getSMILES2(genff,true);
		cout << "Ended getSMILES\n";
		logfile << "\n\n"<<RCOUNT<< ":\tReceived SMILES: "<<get<0>(res) << "\n";
		std::string synthdata;
		#ifndef NOZMQ
			cout <<"ZMQ Code\n";
			//zmqio::ZMQIOThread mythread(logfile);
			//zmq::socket_t& mysock=mythread.getSocket();
			zmq::socket_t& mysock=zmqio::getSocket(src);
			//zmqio::init(logfile);
			//zmq::socket_t& mysock=zmqio::getZMQSocket();
			std::string new_charges=zmqio::communicateWithSocket(mysock,get<0>(res),logfile);
			if(SYNTHGET && new_charges.find("Bad request! check valence")==std::string::npos) synthdata=zmqio::communicateWithSocket(mysock,"!"+get<0>(res),logfile);
		#else
			netio::init(logfile);
			int mysock = netio::getSocket();
			std::string new_charges = netio::communicateWithSocket(mysock,get<0>(res),logfile);
			if(SYNTHGET && new_charges.find("Bad request! check valence")==std::string::npos) synthdata=netio::communicateWithSocket(mysock,"!"+get<0>(res),logfile);
		#endif
		logfile.flush();
		if(new_charges.find("Bad request! check valence")!=std::string::npos)
		{
			if(!FFCHARGES)
			{
				badsmileslog << get<0>(res) << "\n";
				testm->dumpMol("badsmiles_dump_"+to_string(RCOUNT++)+".pdb",genff); badsmileslog.flush();
			}
			delete testm; failures[src]++; testm=nullptr;
			return;
		}
		else
		{
			/*auto splt=stringfx::split(new_charges,',');
			for(std::string& s : splt) cout << s << "\t"; cout << "\n";*/
			try{molfx::assignChargesFromData(txtio::loadChargesFromJSONString(new_charges),get<1>(res),testm,&(genff->getCategory("halogen")));}
			catch(UnmatchingSMILESAndStructurePairException& ex)
			{
				cout << "Bad smiles. Data length is shorter than atom-count\n";
				if(!FFCHARGES)
				{
					badsmileslog << get<0>(res) << "\n";
					testm->dumpMol("badsmiles_dump_"+to_string(RCOUNT++)+".pdb",genff);
					badsmileslog.flush();
				}
				delete testm; failures[src]++; testm=nullptr;
				return;
			}
			testm->dumpMolAsTxt("with_charges_"+to_string(netcomm::RCOUNT++)+".txt");
			if(envprot) envprot->calculateNonBondingEnergy(testm);
			if(SYNTHGET)
			{
				if(SMOOTHSYNTHSCORES) testm->pE=SYNTHPENALTY*(1-std::stod(synthdata));
				else
				{
					int cls=std::stoi(synthdata);
					if(!cls) testm->pE=SYNTHPENALTY;
				}
			}
		}
		cout << "Request processed successfully\n";
	}

}
//Completing the forces namespace
namespace forcefx
{
	static inline Eigen::Vector3d getLengthForce(Atom* a1,Atom* a2,const std::pair<double,double>& lp)
	{
		Eigen::Vector3d r=a1->bondVectorTo(a2);
		if(::abs(r.norm()-get<0>(lp))<0.01) return Eigen::Vector3d(0,0,0);
		return get<1>(lp)*(get<0>(lp)-r.norm())*(r/r.norm());
	}
	static inline Eigen::Vector3d getAngleForce(Atom* a1,Atom* a2,Atom* a3,const std::vector<double>& ap)
	{
		double a=chemtools::getAngle(a1,a2,a3);
		//cout << "Angle: "<<a<<" for equil "<<ap[0]<<"; TOL = "<<(ANGVAR/2.5)<<"\n";
		//if(::abs(a-ap[0])<ANGVAR/5) return Eigen::Vector3d(0,0,0);
		Eigen::Vector3d r=a1->bondVectorTo(a3),s=a2->bondVectorTo(a3); double sn=s.norm();
		Eigen::Vector3d f=r-(r.dot(s)/(sn*sn))*s; f/=f.norm();
		//cout << "Angmag: "<<Eigen::Vector3d((ap[1]/sn)*(ap[0]-a)*f+ap[3]*(ap[2]-r.norm())*(r/r.norm())).norm()<<"\n";
		return (ap[1]/sn)*(ap[0]-a)*f+ap[3]*(ap[2]-r.norm())*(r/r.norm());
	}
	static inline Eigen::Vector3d getDihedralForce(Atom* a1,Atom* a2,Atom* a3,Atom* a4,const std::vector<double>& dp)
	{
		double a=chemtools::getDihedral(a1,a2,a3,a4);
		//cout << "Dihedral: "<<a<<" for equil "<<dp[0]<<"with multiplicity: "<<dp[2]<<"\n";
		double minD=dp[0],ST=(2*PI)/dp[2]; while(minD>=ST-0.05) minD-=ST;
		double d=sin(dp[2]*a-minD); //if(d<1e-1 && cos(dp[2]*a-dp[0])<-0.1) d=throwarandompoint(0,1)<0.5?-1:1;
		Eigen::Vector3d r=a3->bondVectorTo(a4),b=a2->bondVectorTo(a3); double bn=b.norm(),rn;
		Eigen::Vector3d f=b.cross(r),c=a1->bondVectorTo(a2).cross(a2->bondVectorTo(a3));
		//if(quickgeom::angleBetween(f,c)<PI/2) d=-d;
		//cout << d << "\n";
		r=r-(r.dot(b)/(bn*bn))*b;
		rn=r.norm();
		f/=f.norm();
		//cout << ((d*dp[1])/rn) << "\n";
		return ((dp[1]*d)/rn)*f;
	}
	static inline std::vector<Eigen::Vector3d> getImproperForce(Atom* a1,Atom* a2,Atom* a3,Atom* a4,const std::pair<double,double>& ip)
	{
		double a=chemtools::getImproper(a1,a2,a3,a4);
		if(::abs(get<0>(ip)-a)<DIHVAR/2.5 || PI-::abs(get<0>(ip)-a)<DIHVAR/2.5)
		{
			std::vector<Eigen::Vector3d> fc(1,Eigen::Vector3d(0,0,0)); fc.push_back(Eigen::Vector3d(0,0,0)); fc.push_back(Eigen::Vector3d(0,0,0));
			return fc;
		}
		double f=1;
		if(a>PI/2) {a=PI-a; f=-1;}
		f*=get<1>(ip)*(get<0>(ip)-a);
		//double f=get<1>(ip)*(get<0>(ip)-a);
		Eigen::Vector3d b12=a1->bondVectorTo(a2),b13=a1->bondVectorTo(a3),b14=a1->bondVectorTo(a4); //,C=(a2->getPosition()+a3->getPosition()+a4->getPosition())/3;
		Eigen::Vector3d p123=b12.cross(b13),p134=b13.cross(b14),p142=b14.cross(b12); //,plane=a2->bondVectorTo(a3).cross(a3->bondVectorTo(a4)),perp=C-a1->getPosition();
		if(quickgeom::angleBetween(a1->bondVectorTo(a4),p123)>PI/2) p123=-p123;
		if(quickgeom::angleBetween(a1->bondVectorTo(a2),p134)>PI/2) p134=-p134;
		if(quickgeom::angleBetween(a1->bondVectorTo(a3),p142)>PI/2) p142=-p142;
		p123*=f/p123.norm(); p134*=f/p134.norm(); p142*=f/p142.norm();
		//if(quickgeom::angleBetween(plane,perp)>PI/2) plane=-plane;
		std::vector<Eigen::Vector3d> ret; ret.push_back((p123+p142)/(2*b12.norm())); ret.push_back((p123+p134)/(2*b13.norm())); ret.push_back((p134+p142)/(2*b14.norm()));
		return ret;
	}
}
//Completing the Topology class constructor
Topology::Topology(Molecule* m,const ForceField& inpff)
{
	ff=&inpff;
	std::vector<Atom*> atoms=m->getAtoms();
	std::vector<std::vector<int>> bonddata=m->getBonds();
	std::vector<Atom*> visited;
	for(int i=0;i<atoms.size();i++)
	{
		visited.push_back(atoms[i]);
		for(int n1 : bonddata[i])
		{
			if(!contains(visited,atoms[n1]))
			{
				bonds.push_back(make_pair(i,n1));
				bondpars.push_back(ff->getLengthParameters(atoms[i]->toString(),atoms[n1]->toString()));
			}
			for(int n2 : bonddata[n1])
			{
				if(n2==i) continue;
				if(!contains(visited,atoms[n2]))
				{
					std::vector<int> ang(1,i); ang.push_back(n1); ang.push_back(n2);
					angles.push_back(ang);
					anglepars.push_back(ff->getAngleParamtersFailsafe(atoms[i],atoms[n1],atoms[n2]));
				}
				for(int n3 : bonddata[n2])
				{
					if(n3==n1) continue;
					if(contains(visited,atoms[n3])) continue;
					std::vector<int> dih(1,i); dih.push_back(n1); dih.push_back(n2); dih.push_back(n3);
					dihedrals.push_back(dih);
					dihpars.push_back(ff->getDihedralParametersFailsafe(atoms[i],atoms[n1],atoms[n2],atoms[n3]));
				}
				for(int n3 : bonddata[n1])
				{
					if(n3==n2 || n3==n1) continue;
					if(contains(visited,atoms[n3])) continue;
					std::vector<int> imp(1,n1); imp.push_back(i); imp.push_back(n2); imp.push_back(n3);
					impropers.push_back(imp);
					improppars.push_back(ff->getImproperParametersFailsafe(atoms[n1],atoms[i],atoms[n2],atoms[n3]));
				}
			}
		}
		this->atoms.push_back(atoms[i]);
		bool found=false; for(Atom* a : atomtypes) {if(a->toString()==atoms[i]->toString()) {found=true;break;}}
		if(!found) atomtypes.push_back(inpff.getAtom(atoms[i]));
	}
}
std::vector<Eigen::Vector3d> Topology::calculateInternalForces(Molecule* m,bool takeimprop) const
{
	std::vector<Atom*> allats=(m)?m->getAtoms():this->atoms;
	std::vector<Eigen::Vector3d> forces; //,poses;
	//Eigen::Vector3d COM=m->getCentreOfMass();
	Eigen::Vector3d F;
	for(int i=0;i<atoms.size();i++) forces.push_back(Eigen::Vector3d(0,0,0));
		//poses.push_back(atoms[i]->getPosition()-COM);
	for(int i=0;i<bonds.size();i++)
	{
		int c1=get<0>(bonds[i]),c2=get<1>(bonds[i]);
		forces[c2]+=forcefx::getLengthForce(atoms[c1],atoms[c2],bondpars[i]);
		forces[c1]+=forcefx::getLengthForce(atoms[c2],atoms[c1],bondpars[i]);
	}
	for(int i=0;i<angles.size();i++)
	{
		int c1=angles[i][0],c2=angles[i][1],c3=angles[i][2];
		F=forcefx::getAngleForce(atoms[c1],atoms[c2],atoms[c3],anglepars[i]);
		forces[c3]+=F; //forces[c2]-=F;
		F=forcefx::getAngleForce(atoms[c3],atoms[c2],atoms[c1],anglepars[i]);
		forces[c1]+=F; //forces[c2]-=F;
	}
	for(int i=0;i<dihedrals.size();i++)
	{
		int c1=dihedrals[i][0],c2=dihedrals[i][1],c3=dihedrals[i][2],c4=dihedrals[i][3];
		F=forcefx::getDihedralForce(atoms[c1],atoms[c2],atoms[c3],atoms[c4],dihpars[i]);
		forces[c4]+=F; //forces[c1]-=F;
		F=forcefx::getDihedralForce(atoms[c4],atoms[c3],atoms[c2],atoms[c1],dihpars[i]);
		forces[c1]+=F; //forces[c4]-=F;
	}
	if(takeimprop)
	{
		for(int i=0;i<impropers.size();i++)
		{
			if(!get<1>(improppars[i])) continue;
			int c1=impropers[i][0],c2=impropers[i][1],c3=impropers[i][2],c4=impropers[i][3];
			std::vector<Eigen::Vector3d> foc=forcefx::getImproperForce(atoms[c1],atoms[c2],atoms[c3],atoms[c4],improppars[i]);
			forces[c2]+=foc[0]/3;
			forces[c3]+=foc[1]/3;
			forces[c4]+=foc[2]/3;
			forces[c1]-=(foc[0]+foc[1]+foc[2])/3;
		}
	}
	std::vector<std::pair<int,int>> ignore=bonds;
	for(int i=0;i<angles.size();i++) ignore.push_back(make_pair(angles[i][0],angles[i][2]));
	for(int i=0;i<dihedrals.size();i++) ignore.push_back(make_pair(dihedrals[i][0],dihedrals[i][3]));
	//Non-bonding part
	/*for(int i=0;i<atoms.size();i++)
	{
		for(int j=i+1;j<atoms.size();j++)
		{
			bool found=false; for(const auto& p : ignore) {if((get<0>(p)==i && get<1>(p)==j) || (get<0>(p)==j && get<1>(p)==i)) {found=true; break;}}
			if(found) continue;
			Eigen::Vector3d r=atoms[j]->bondVectorFrom(atoms[i]); double n=r.norm();
			if(n<0.25) continue;
			//cout << n << "\n";
			r=(chemtools::getLJForce(atoms[i],atoms[j],n)/n)*r+(chemtools::getElectrostaticForce(atoms[i],atoms[j],n)/n)*r;
			forces[i]+=r;
			forces[j]-=r;
		}
	}*/
	return forces;
}
double Topology::calculateTotalAngleEnergy(Molecule* src) const
{
	double E=0;
	std::vector<Atom*> ats=(src)?src->getAtoms():this->atoms;
	for(int i=0;i<anglepars.size();i++) E+=angleEnergyFunction(anglepars[i],chemtools::getAngle(ats[angles[i][0]],ats[angles[i][1]],ats[angles[i][2]]),ats[angles[i][0]]->distanceFrom(ats[angles[i][2]]));
	return E;
}
double Topology::calculateTotalDihedralEnergy(Molecule* src) const
{
	double E=0;
	std::vector<Atom*> ats=(src)?src->getAtoms():this->atoms;
	for(int i=0;i<dihpars.size();i++) E+=dihedralEnergyFunction(dihpars[i],chemtools::getDihedral(ats[dihedrals[i][0]],ats[dihedrals[i][1]],ats[dihedrals[i][2]],ats[dihedrals[i][2]]));
	return E;
}

//Completing physics (namespace)
namespace physics
{
	void gradientDescentPhysics(const System& sys,Molecule* target,int nsteps,double D,double R,double L,double L2)
	{
    D/=target->getMass();
		double MOI=target->getMomentOfInertia();
		cout << "System has moment of inertia: "<<MOI<<"\n";
		if(MOI!=MOI) return;
		R/=MOI;
    std::vector<std::pair<std::vector<Eigen::Vector3d>,std::vector<Eigen::Vector3d>>> forces;
		double totbe=0;
    for(int stepcount=0;stepcount<nsteps;stepcount++)
    {
			target->myBE=0;
      for(Molecule* am : sys.getMolecules())
      {
        if(am==target) continue;
        forces.push_back(am->calculateActingForces(target));
		    for(Atom* ta : target->getAtoms()) target->myBE+=get<0>(ta->setEnergyContribution(am->calculateNonBondingEnergy(ta)));
      }
      if(!forces.size()) {cerr << "WARN: No forces found! Are there no molecules other than target in the system for gradient descent?\n"; continue;}
      std::vector<Eigen::Vector3d> sumf=get<0>(forces[0]),atposes=get<1>(forces[0]);
      for(int i=1;i<forces.size();i++) {for(int j=0;j<sumf.size();i++) sumf[j]+=get<0>(forces[i])[j];}
      std::pair<std::vector<Eigen::Vector3d>,std::vector<Eigen::Vector3d>> totalforces=make_pair(sumf,atposes);
			Eigen::Vector3d COM=target->getCentreOfMass();

      Eigen::Vector3d totf(0,0,0); for(const Eigen::Vector3d& v : get<0>(totalforces)) totf+=v;
			Eigen::Vector3d tott(0,0,0); for(const Eigen::Vector3d& v : get<1>(totalforces)) tott+=v;
      Eigen::Vector3d shift=totf*D; double mdist=shift.norm();
			Eigen::Vector3d rot=tott*R; double rn=rot.norm();
			if(mdist<L && rn<L2)
			{
				if(!QUIET) cout << "Step size small enough. Will stop gradient descent.\n";
				return;
			}
      if(!QUIET) cout << "Step: "+to_string(stepcount)+" with force "+to_string(totf.norm())+"*1e-9 kN/mol with shift: "+to_string(mdist)+" with potential energy: "+to_string(target->myBE)+" kJ/mol. Acting torque: "+to_string(tott.norm())+" with rot: "+to_string(rn)+"\n";
      for(Atom* a : target->getAtoms()) a->setPosition(a->getPosition()+shift);
			if(!QUIET) cout << "Moved!";
				//cout << shift.transpose() << "\n";
			target->rotateMolecule(rot,COM); //Magnitude of vector rot is used to calculate angle. The normal vector direction is the unit vector in this direction
			if(!QUIET) cout << "Rotated!";
				//cout << rot.transpose() << "\n";
			forces=std::vector<std::pair<std::vector<Eigen::Vector3d>,std::vector<Eigen::Vector3d>>>();
    }
		if(!QUIET) cout << "Total number of steps required were completed. Stopping gradient descent.\n";
  }
	void atomisticGradientDescentPhysics(const System& sys,Molecule* target,int nsteps,double D,double DU,double DS,double R,double L)
	{
		Topology mytop(target,*sys.ff);
		std::vector<Atom*> atoms=target->getAtoms();
		std::vector<Eigen::Vector3d> forces; for(int i=0;i<target->getSize();i++) forces.push_back(Eigen::Vector3d(0,0,0));
		double molmass=target->getMass();
		for(int stepcount=0;stepcount<nsteps;stepcount++)
    {
			//if(stepcount==(nsteps*1)/5) target->dumpMol("inframe.pdb",sys.ff);
			for(Molecule* am : sys.getMolecules())
      {
				std::vector<Eigen::Vector3d> tf;
        if(am==target)
				{
					//if(stepcount>0.2*nsteps) { tf=mytop.calculateInternalForces(target,false); cout << "No impropers in step "<<stepcount<<"\n";}
					tf=mytop.calculateInternalForces(target,true);
				}
				else tf=get<0>(am->calculateActingForces(target));
				for(int i=0;i<tf.size();i++) forces[i]+=tf[i]*((am==target)?D:R);
			}
			if(!forces.size()) {cerr << "WARN: No forces found! Are there no molecules other than target in the system for gradient descent?\n"; continue;}
			double fM=0; int fI=-1;
			for(int i=0;i<forces.size();i++)
			{
				if(fM<forces[i].norm())
				{
					fM=forces[i].norm();
					fI=i;
				}
			}
			if(fI==-1) {cout << "Warn: Forces hit zero... Check the molecule for confirmation\n"; return;}
			if(!QUIET) cout << "Chose atom at index: "<<fI<<" for maximum force value: "<<fM<<" *1e-9 kN/mol\n";
			if(atoms[fI]->isHydrogen()) forces[fI]/=12; else forces[fI]/=atoms[fI]->seek_mass();
			double mdist=forces[fI].norm();
			//double mdist=fM/atoms[fI]->seek_mass();
			if(mdist<L)
			{
				if(!QUIET) cout << "Step size small enough("<<L<<"). ";
				if(D<DU) D+=DS;
				else
				{
					if(!QUIET) cout <<" Will stop gradient descent.\n";
					return;
				}
			}
			//else if(mdist>0.0015) {mdist=0.0015; forces[fI]*=(mdist/forces[fI].norm());}
			if(!QUIET) cout << "Step: "+to_string(stepcount)+" with force "+to_string(fM)+"*1e-9 kN/mol with shift: "+to_string(mdist)+"\n";
			Eigen::Vector3d oldpos=atoms[fI]->getPosition();

			if(forces[fI].norm()>0.05)
			{
				forces[fI]*=atoms[fI]->seek_mass()/molmass;
				for(Atom* ia : target->getAtoms()) ia->setPosition(ia->getPosition()+forces[fI]);
			}
			else atoms[fI]->setPosition(atoms[fI]->getPosition()+forces[fI]);

			if(!QUIET) cout << "Moved!";
				//if(stepcount%250==0) target->dumpMol("inprogress_mol_"+to_string(stepcount)+".pdb",sys.ff);
		}
	}
	void optimizeMolecule(Molecule* m,const ForceField& ff,int nsteps,double D)
	{
		System temp(ff); temp.addMolecule(m);
		atomisticGradientDescentPhysics(temp,m,nsteps,D);
	}
}
#endif
