#ifndef HPP_ATOM
#define HPP_ATOM 1
#include<iostream>
#include<cstring>
#include "support/externio.hpp"

using namespace std;

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#########################Define the ATOM class##############################|
|##################To store atom/particle parameters#########################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/

//Exceptions
/**\deprecated This exception will be thrown if any atom object does not have it's "atom_type" uninitialized*/
class EmptyAtomTypeException: public exception {}; //When Atom type variable is not initialized in an Atom object
/**This exception is thrown if any attempt is made to add a bond to an atom with no free valencies. See Atom::seal() and Atom::free() for managing valencies*/
class NoValencyAvailableException: public exception {}; //Self-explanatory
class Molecule;
class MolecularFragment;
class ForceField;
class System;

typedef std::vector<MolecularFragment*> FragmentSet;
typedef std::vector<Eigen::Vector3d> Pose;

/** \class Atom Atom.hpp "graph/Atom.hpp"
	This is the Atom class. It is the most fundamental building-block for molecule generation using DeNovo.
	An atom object contains all the force-field parameters (see: ForceField) and position data contained in it.
	The Molecule class uses arrays of these Atom* (pointers) and defines connections between them in order to manage molecules
	All position information is stored in double precision and the unit of distance in DeNovo is nm

	This class makes use of the Eigen libraries. See Eigen: http://eigen.tuxfamily.org/ for more details. (git repo: https://gitlab.com/libeigen/eigen.git)
*/
class Atom
{
	Molecule* myMol=nullptr;
	int valency=0,hybridization=3,standval,resnum=1;
	std::string atom_type="DUM",residue="DNV";
	//char* atom_type;
	double x,y,z,charge,sigma,epsilon,vanderwaal_radius,mass,enercont=0,hsenercont=0;
	char chain='A';
	bool staticcycl=false;
	double formalcharge=0;
	float solvwt=0;
	public:
		bool seed=false; /**< Is this a seed atom? */
		bool cycl=false; /**< Has the atom cyclized? */
		bool mustCycl=false; /**< Does this atom need to cyclize? */
		int hashvalue=0; /**<\deprecated Atom's unique hash value (deprecated - See Atom::hashvector) */
		std::vector<int> hashvector; /**< Atom's unique hash-vector (still under development) */
		double posdev=0; /**<Deviation from equilibrium position - Used to calculate error in potential energy calculations (Default: 0)*/
		bool genpurpflag=false; /**< General purpose flag */
		int bindpoints=0;
		#ifdef USEGENERALDATA
		int genpurpint=0; /**<General purpose integer storage*/
		#endif
	public:
		/**@brief Blank default constructor - Try not to use this */
		Atom();
		/**@brief Construct an atom with only position and valency
			 @details Construct a dummy "DUM" atom at a given location. It is given valency 'v' (default: 0). It carries L-J parameters close to a Carbon atom.<br/> It is usually used as a placeholder*/
		Atom(double x,double y,double z,int v=0) : Atom("DUM",x,y,z,v) {}
		/**@brief Construct an atom with only name, position and valency
			 @details Construct an atom with given name at a given location. It is given valency 'v' (default: 0). It carries L-J parameters close to a Carbon atom.<br/> It is usually used as a placeholder*/
		Atom(const std::string& name,double x,double y,double z,int v=0);
		/**@brief The main Atom constructor.
			 @details This constructor loads all the data required including charge, L-J parameters (sigma,epsilon) as well as position (x,y,z) and valency (v).<br/> It is recommended you use this constructor if you ever have to directly construct an atom. See also: ForceField::getAtom() */
		Atom(const std::string& name,double x,double y,double z,int v,double weight,double charge,double sigma,double epsilon,double vanderwaal_rad);
		/** @brief Another similar constructor. This also allows specifying the hybridization
			  @details Another main Atom constructor. See Atom(char* name,double x,double y,double z,int v,double weight,double charge,double sigma,double,double) .<br/> This constructor also allows you to explicitly mention the hybridization of this atom. See also: ForceField::getAtom()*/
		Atom(const std::string&,double,double,double,int,double,double,double,double,double,int);
		/**@brief The standard copy constructor.
			 @details The standard copy constructor is safe to use directly (barring some technicalities).<br/> Note the definition. You cannot pass a pointer to copy. Please dereference it ONLY before passing to this method */
		Atom(const Atom&);
		/**@brief The standard destructor for the Atom class.*/
		~Atom();
		void display_atomtype();
		void display_position();
		std::string seek_atomtype() const;
		/**@brief Get the x-coordinate from the atom's position*/
		inline double seek_x() const;
		/**@brief Get the y-coordinate from the atom's position*/
		inline double seek_y() const;
		/**@brief Get the z-coordinate from the atom's position*/
		inline double seek_z() const;
		/**@brief Set the x-coordinate of the atom's position*/
		inline void modify_x(double);
		/**@brief Set the y-coordinate of the atom's position*/
		inline void modify_y(double);
		/**@brief Set the z-coordinate of the atom's position*/
		inline void modify_z(double);
		/**@brief Get the number of free valencies this atom has */
		inline int seek_valency() const;
		/**@brief Get the force-field (or otherwise loaded) charge parameter for this atom*/
		inline double seek_charge() const;
		/**@brief Get the force-field (or otherwise loaded) sigma parameter for this atom*/
		inline double seek_sigma() const ;
		/**@brief Get the force-field (or otherwise loaded) epsilon parameter for this atom*/
		inline double seek_epsilon() const ;
		/**@brief Get the force-field (or otherwise loaded) vanderwaal radius for this atom*/
		inline double seek_radius() const;
		/**@brief Get the mass of this atom*/
		inline double seek_mass() const {return mass;}
		/**@brief Gives the atom name for this atom*/
		inline const std::string& toString() const {return atom_type;}
		/**@brief Make the atom a seed
		 @details A wrapper for "seed=true". This atom is then declared to be a seed atom. See the algorthim for more details*/
		inline void makeSeed() {seed=true;}
		/**@brief A wrapper for "seed". Checks if this atom is a seed or not*/
		inline bool isSeed() {return seed;}
		/**@brief Set the number of free valencies
		@details Allows manually setting the number of free valecies for this atom. Unless you are using an unparametrized atom-type (or are not using a Force-field at all), it is not recommended to use this method*/
    inline void setValency(int n) {valency=n;}
		/**@brief Get the force-field (or otherwise loaded) maximum allowed valency for this atom*/
    inline int getStandardValency() const {return standval;}
		/**@brief Get the approximate solvation correction for this atom-type. It must be preloaded from the ForceField*/
		inline float getSolvationCorrection() const {return solvwt;}
		/**@brief Set the approximate solvation correction for this atom-type.
			 @details Use to manually set the contribution weight for this atom-type. If you are doing this manually for one atom-type, you are likely going to have to do it for all in the ForceField*/
		inline void setSolvationCorrection(float w) {solvwt=w;}

		//Methods added by Sreyas (13-12-2019 onwards)
		/**@brief Returns this complete position of the enitre atom in an "Eigen::Vector3d" object.*/
		inline Eigen::Vector3d getPosition() const;
		/**@brief Allows setting the position of the atom by providing an "Eigen::Vector3d" object.*/
		inline void setPosition(const Eigen::Vector3d& v);
		/**@brief Allows setting the position of the atom by providing three coordinates (x,y,z).*/
		inline void setPosition(const double& x,const double& y, const double& z);
		/**@brief Get a vector from this atom <b>to</b> 'a'
		  @details Returns an Eigen::Vector3d object pointing from "this" atom <b>to the target atom</b> (See also: bondVectorFrom(Atom*) const).<br/> The magnitude of this vector will match the value of distanceFrom(Atom* a)
			@return Vector (Eigen::Vector3d object)
		*/
		inline Eigen::Vector3d bondVectorTo(Atom* a) const;
		/**@brief Get a vector <b>from</b> 'a' to this atom
			 @details Returns an Eigen::Vector3d object pointing to "this" atom <b>from the target atom</b> (See also: bondVectorTo(Atom*) const).<br/>The magnitude of this vector will match the value of distanceFrom(Atom* a)
			 @return Vector (Eigen::Vector3d object)
		*/
		inline Eigen::Vector3d bondVectorFrom(Atom* a) const;
		/**@brief Evaluates the distance from "this" atom to another
			 @return double: distance (in nm)*/
		inline double distanceFrom(const Atom* a2) const;
		/**@brief Evaluates the distance from "this" atom to a random point specified by the vector 'v'
			 @return double: distance (in nm)
		*/
		inline double distanceFrom(const Eigen::Vector3d& v) const;
		/**@brief Tests if this atom is a hydrogen atom
			 @details Hydrogen atoms are defined to be'light' atoms. This method is used to make the distinction between light and heavy atoms. */
		inline bool isHydrogen() const {return atom_type[0]=='H';}
		/**@brief Get the force-field (or otherwise loaded) hybridization of this atom*/
		inline int getHybridization() const {return hybridization;}
		/**@brief A simple wrapper to reduce the valency by 1*/
		inline void addBond() {valency--;}
		/**@brief Checks if this atom can make further bonds
			 @details This method basically checks if the atom has any free valencies or not. It is equivalent to seek_valency()>0*/
		inline bool canBond() {return valency;}
		/**@brief Store the source residue name in this Atom object
			@details Sets the residue name for this Atom object to that of the residue it belongs to.<br/> It is useful when loading proteins, dna, etc. from a forcefield with specific atom names for atom of each of it's residues.<br/> By default no standard residue name is taken (default residue name is DNV).*/
		inline void setResidue(const std::string& r) {residue=r;}
		/**@brief Store the source residue number in this Atom object
			@details Sets the residue number for this Atom object to that of the residue it belongs to.<br/> It is useful when loading proteins, dna, etc. from a forcefield with specific atom names for atom of each of it's residues.<br/> (Default residue number is 1).*/
		inline void setResidueNumber(const int r) {resnum=r;}
		/**@brief Get the assigned residue name
			 @details Get back the residue name which was set for this atom.<br/>
			 Residue names are automatically loaded when loading from a PDB or GRO file (see  Molecule::Molecule(const std::string&,const ForceField&,const std::string&,bool,bool,char))<br/>Alternatively, then can be set by: setResidue(const std::string& s)
			 @return const std::string&: The residue name for this atom
		*/
		inline const std::string& getResidue() const {return residue;}
		/**@brief Get the assigned residue number
			 @details Get back the residue name which was set for this atom.<br/>
			 Residue numbers are automatically loaded when loading from a PDB or GRO file (see  Molecule::Molecule(const std::string&,const ForceField&,const std::string&,bool,bool,char))<br/>Alternatively, then can be set by: setResidue(const std::string& s)
			 @return const std::string&: The residue name for this atom
		*/
		inline int getResidueNumber() const {return resnum;}
		/**@brief Returns the Molecule* (pointer) to which it belongs (if assigned)
			 @details If assigned, returns the Molecule to which this atom belongs. Unless you are very sure that you have assigned the pointers of the parent Molecule objects to the atoms, do NOT use this method.<br/>Consider Molecule::getAtom(int), and Molecule::getAtoms().<br/>It is better to access atom pointrs from the Molecule object than vice-versa
			 @return Molecule*: Pointer to the Molecule object in which it belongs
		*/
		inline Molecule* getMolecule() const {return myMol;}
		/**@brief Assign a Molecule to this atom
			 @details Allows you to assign a Molecule (pointer) to this atom. It is usually meant for the tracability of parent Molecule object from the atoms.<br/>It is <b>not recommended</b> that you assign any other Molecule to this atom as it might not be compatible with many standard procedures in DeNovo*/
		inline void setMolecule(Molecule* mp) {myMol=mp;}
		/**@brief Do NOT use this method directly. It manually sets the interaction energy of this atom to the Protein (or other macromolecule)
			 @details This method may cause severe errors in computing "DeNovo interaction energy" for growth. See: Protein::calculateNonBondingEnergy(Molecule*).<br/> This method assigns each atom's energy contribution*/
		inline double setEnergyContribution(double e) {return enercont=e;}
		/**@brief Set energy contribution value of total and hotspot together as a pair (total,hotspot). See Also: setEnergyContribution(double), setHotspotEnergyContribution(double)*/
		inline std::pair<double,double> setEnergyContribution(const std::pair<double,double>& p) {enercont=get<0>(p); hsenercont=get<1>(p); return p;}
		/**@brief Get the interaction energy of this atom with the Protein (or other macromolecule in the System)
			 @details Get the calculated energy contribution of this atom. If not calculated (by Protein::calculateNonBondingEnergy(Molecule *) or otherwise), it is 0*/
		inline double getEnergyContribution() const {return enercont;}
		/**@brief Do NOT use this method directly. It manually sets the interaction energy of this atom to the "hotspot"
			 @details This method may cause severe errors in computing "DeNovo interaction energy" for growth. See: Protein::calculateNonBondingEnergy(Molecule*).<br/> This method assigns each atom's energy contribution (even to the hotspot). See Protein class for explanation of "hotspot"s*/
		inline double setHotspotEnergyContribution(double e) {return hsenercont=e;}
		/**@brief Get the interaction energy of this atom to the hotspot.
			 @details Get the calculated energy contribution of this atom to the hotspot (see Protein class) atoms. If not calculated (by Protein::calculateNonBondingEnergy(Molecule*) or otherwise), it is 0*/
		inline double getHotspotEnergyContribution() const {return hsenercont;}
		/**@brief Get the chain ID (as a character) if specified.
			 @details This method is used to get the Chain ID to which this atom belongs. It is useful when using Proteins with multiple chains.<br/>
			 Default: 'A'
			 @return The one-character ID of the chain (as according to PDB)*/
		inline char getChain() const {return chain;}		
		/**@brief Specify the original chain ID from which this atom was picked.
			 @details: Allows user to manually specify the chain ID (in case multiple chains are present). If nothing is specified, chain is 'A'<br/>See also: Molecule::Molecule(const std::string&,const ForceField&,const std::string&,bool,bool,char)*/
		inline void setChain(char c) {chain=c;}
		/**@brief A wrapper for cycl=true/false*/
		inline bool setCyclized(bool v=true) {return cycl=v;}
		/**@brief A wrapper for cylc. Checks is the atom has cyclized or not*/
    inline bool isCyclized() const {return cycl;}
		/**@brief Do not allow this atom to make further bonds
			@details This method rejects all remaining valencies of this atom. The number of free valecies of this atom are set to 0 (irrespective of it connection count and max-valency).<br/>It is useful when using polyatomic seeds (see algorithm \ref test.html "PDF" for details about seeds) if you want to grow from a specific atom only*/
		inline void seal() {valency=0;}
		/**@brief Manually set partial charge for this Atom object
			 @details Allows rewriting the initally loaded charge parameter for the atom. It is highly recommended that you load the atom with proper parameters using complete force-fields (see: ForceField), and use ForceField::getAtom(const std::string& atom_name) to get the complete atom (with all information)<br/>
			  Do not load the atom first and then assign the charge. Also, manual charge assignment is rarely a good substitute for well parameterized force-fields when using DeNovo*/
		inline void setCharge(double c) {charge=c;}
		/**@brief Developer only method. Allows "renaming" an atom.
			 @details This method will change the "name" of the atom. Note that none of the parameters are changed even if the new name is defined differently in the force-field.*/
		inline void rename(const std::string& an) {atom_type=an;}
		/**@brief Set the (simulated or otherwise obtained) deviation from mean position of the atom*/
		inline void setPositionDeviation(double v) {posdev=v;}
		/**@brief Get the (simulated or otherwise obtained) deviation from mean position of the atom*/
		inline double getPositionDeviation() const {return posdev;}
		/**@brief Preserve cyclization state (cyclized? y/n) when copying this atom (Default: off)
			 @details This can be useful if you want to set the cyclization state of a template molecule's atoms before generation*/
		inline void makeCycleStatic() {staticcycl=true;}
		/**@brief Is cyclization state (cyclized? y/n) preserved when copying this atom (Default: off). See Also: makeCycleStatic()*/
		inline bool isCycleStatic() const {return staticcycl;}
		/**@brief Assign the standard valency value for this atom*/
		inline void setStandardValency(int sv) {standval=sv;}

		 /**@brief Return the (loaded) formal charge for this atom type (this is usually 0)*/
		 inline double getFormalCharge() const {return formalcharge;}
		 /**@brief Assign (load) a formal charge for this atom*/
		 inline void setFormalCharge(double fc) {formalcharge=fc;}
		/**@brief Reassign basic parameters (including atom-type name) to match supplied atom*/
		inline void assignAtom(const Atom* to_copy,bool copypos=false)
		{
		  standval=to_copy->standval;
			solvwt=to_copy->solvwt;
			//if(atom_type) delete[] atom_type;
			atom_type=to_copy->atom_type;
			if(copypos)
			{
				x=to_copy->x;
				y=to_copy->y;
				z=to_copy->z;
			}
			charge=to_copy->charge;
			sigma=to_copy->sigma;
			epsilon=to_copy->epsilon;
			vanderwaal_radius=to_copy->vanderwaal_radius;
			mass=to_copy->mass;
			hybridization=to_copy->hybridization;
			mustCycl=to_copy->mustCycl;
			formalcharge=to_copy->formalcharge;
		}
	private:
		void evaluateHybridization();
};

/*
Define class member functions here
*/

Atom::Atom(){
	valency=4;
  standval=valency;
	//atom_type=new char[5];
	atom_type = "DUM";
	//strcpy(atom_type,typeat.c_str());
	x=0.0;y=0.0;z=0.0;
	mass=12.065;
	charge=0.0;
	sigma=0.405358916754; epsilon = 0.08368;
	vanderwaal_radius = sigma/2.0;
	cycl=false;
}

Atom::Atom(const std::string& nm, double xx, double yy, double zz, int v){
	valency=v;
  standval=valency;
	atom_type = nm;
	x=xx;y=yy;z=zz;
	mass=12.065;
	charge=0.0;
	sigma=0.405358916754; epsilon = 0.08368;
	vanderwaal_radius = sigma/2.0;
	cycl=false;
}

Atom::Atom(const std::string& id,double xx,double yy,double zz,int v,double wt,double c,double s,double e,double r)
{
	//atom_type=id;
	//atom_type=new char[10];
	atom_type=id;
	//strcpy(atom_type,id);
	x=xx;
	y=yy;
	z=zz;
	valency=v;
  standval=valency;
	mass=wt;
	charge=c;
	sigma=s;
	epsilon=e;
	vanderwaal_radius=r;
	evaluateHybridization();
	cycl=false;
}
Atom::Atom(const std::string& id,double xx,double yy,double zz,int v,double wt,double c,double s,double e,double r,int h) : Atom::Atom(id,xx,yy,zz,v,wt,c,s,e,r) {hybridization=h;}
/*{
	//atom_type=id;
	//atom_type=new char[10];
	//strcpy(atom_type,id);
	atom_type=std::string(id);
	x=xx;
	y=yy;
	z=zz;
	valency=v;
  standval=valency;
	mass=wt;
	charge=c;
	sigma=s;
	epsilon=e;
	vanderwaal_radius=r;
	hybridization=h;
	cycl=false;
}*/
Atom::Atom(const Atom& to_copy)
{
	valency=to_copy.valency;
  standval=to_copy.standval;
	residue=to_copy.residue;
	//if(atom_type) delete[] atom_type;
	atom_type=to_copy.atom_type;
	x=to_copy.x;
	y=to_copy.y;
	z=to_copy.z;
	charge=to_copy.charge;
	sigma=to_copy.sigma;
	epsilon=to_copy.epsilon;
	vanderwaal_radius=to_copy.vanderwaal_radius;
	mass=to_copy.mass;
	hybridization=to_copy.hybridization;
	cycl=false;
	enercont=to_copy.enercont;
	hsenercont=to_copy.hsenercont;
	mustCycl=to_copy.mustCycl;
	staticcycl=to_copy.staticcycl;
	if(to_copy.isCycleStatic()) cycl=to_copy.cycl;
	formalcharge=to_copy.formalcharge;
	solvwt=to_copy.solvwt;
}
Atom::~Atom()
{
	//delete[] atom_type;
}
void Atom::display_atomtype()
{
	cout<<"Atom Type = "<<atom_type<<"\t";
}
void Atom::display_position()
{
	cout<<"X = "<<x<<"\tY = "<<y<<"\tZ = "<<z<<"\tValency = "<<valency<<endl;
}
inline std::string Atom::seek_atomtype() const {return atom_type;}
inline double Atom::seek_x() const {return x;}
inline double Atom::seek_y() const {return y;}
inline double Atom::seek_z() const {return z;}

inline void Atom::modify_x(double x1) {x=x1;}
inline void Atom::modify_y(double y1) {y=y1;}
inline void Atom::modify_z(double z1) {z=z1;}
inline int Atom::seek_valency() const {return valency;}
inline double Atom::seek_radius() const {return vanderwaal_radius;}
inline double Atom::seek_charge() const {return charge;}
inline double Atom::seek_sigma() const {return sigma;}
inline double Atom::seek_epsilon() const {return epsilon;}

//Method bodies added by Sreyas
inline Eigen::Vector3d Atom::getPosition() const {return Eigen::Vector3d(x,y,z);}
inline void Atom::setPosition(const Eigen::Vector3d& v) {x=v(0); y=v(1); z=v(2);}
inline void Atom::setPosition(const double& xv,const double& yv, const double& zv) {x=xv;y=yv;z=zv;}
inline double Atom::distanceFrom(const Atom* a2) const {return sqrt(sqr(seek_x()-a2->seek_x()) + sqr(seek_y()-a2->seek_y()) + sqr(seek_z()-a2->seek_z()));}
inline double Atom::distanceFrom(const Eigen::Vector3d& v) const {return sqrt(sqr(seek_x()-v(0))+sqr(seek_y()-v(1))+sqr(seek_z()-v(2)));}
inline Eigen::Vector3d Atom::bondVectorTo(Atom* a) const {return Eigen::Vector3d(a->seek_x()-seek_x(),a->seek_y()-seek_y(),a->seek_z()-seek_z());}
inline Eigen::Vector3d Atom::bondVectorFrom(Atom* a) const {return Eigen::Vector3d(seek_x()-a->seek_x(),seek_y()-a->seek_y(),seek_z()-a->seek_z());}

void Atom::evaluateHybridization()
{
	//if(!atom_type) {throw EmptyAtomTypeException();}
	switch(atom_type[0])
	{
		case 'C': hybridization=valency-1; break; //4 valency -> sp[3]
		case 'S': case 'O': hybridization=valency+1; break; //1 valency -> sp[2]
		case 'P': case 'N': hybridization=(valency==4)?3:valency; break; //3 valency -> sp[3]
		case 'F': case 'H': hybridization=0; break;
	}
}

std::ostream& operator<<(std::ostream& os,const Atom& a)
{
	os << a.toString() << "\t" << a.seek_x()*10 << "\t" << a.seek_y()*10 << "\t" << a.seek_z()*10 << "\t" << a.getStandardValency() << "\t" << a.seek_mass()<< "\t"<< a.seek_charge()<<"\t"<<a.seek_sigma()<<"\t"<<a.seek_epsilon()<<"\t"<<a.seek_radius();
	return os;
}

class Molecule;
namespace chemstructs
{
	/** \struct Bond Atom.hpp "graph/Atom.hpp"
		 This is a simple struct to store the components of a bond (three atoms).
		 It also calculates the vector of the bond when the component atoms are provided<br/>
		 This class makes use of the Eigen libraries. See Eigen: http://eigen.tuxfamily.org/ for more details. (git repo: https://gitlab.com/libeigen/eigen.git)
	*/
	struct Bond
	{
		Atom *a1; /**<The first atom of the bond*/
		Atom *a2; /**<The second atom of the bond*/
		int order=1; /**<The bond order. This can be specified by the user later as well*/
		Eigen::Vector3d vect; /**<The vector in the direction of this bond with magnitude equal to bond-length (in nm)*/

		/**@brief The simple Bond constructor which uses two bonds to generate this Bond object*/
		Bond(Atom* at1,Atom* at2) {a1=at1; a2=at2; vect=at1->bondVectorTo(at2);}

		/**@brief Get the length of the bond*/
		inline double getLength() const {return vect.norm();}
	};
	/** \struct Angle Atom.hpp "graph/Atom.hpp"
		 This is a simple struct to store the components of an angle (three atoms). The actual order in which the atoms are provided is of no significance as long as the angle subtended is accurately described (eg. The middle atom must be the second atom)<br/>
		 It also calculates the plane of the three atoms and the angle between them when they are provided
	*/
	struct Angle
	{
		Atom *a1; /**<The first atom forming the angle*/
		Atom *a2; /**<The second atom forming the angle*/
		Atom *a3; /**<The third atom forming the angle*/
		double angle=0; /**<The angle calculated between the component atoms*/
		Eigen::Vector3d plane; /**<The normal vector to the plane of this angle*/

		/**@brief The simple Angle constructor which uses three atoms to construct the angle*/
		Angle(Atom* at1,Atom* at2,Atom* at3)
		{
			a1=at1; a2=at2; a3=at3;
			angle=quickgeom::angleBetween(at2->bondVectorTo(at1),at2->bondVectorTo(at3));
			plane=quickgeom::planeOf(at2->bondVectorTo(at1),at2->bondVectorTo(at3));
		}
	};
}
/**The netcomm nammespace*/
namespace netcomm
{
	static const ForceField* genff=nullptr;
	int RCOUNT=0;
	std::ofstream badsmileslog;
	class MolRequest;
	static int* failures;
	inline void addRequest(System m,int id);// {molq.push(MolRequest(m,id));}
	inline void processRequest(std::ostream&,Molecule*&,int src=0);
	inline int queueSize();
	inline std::vector<MolRequest> drainQueue();
}
namespace chemtools
{
	/**@brief Calculate the angle between any three atoms*/
	static inline double getAngle(Atom* a1,Atom* a2,Atom* a3) {return quickgeom::angleBetween(a2->bondVectorTo(a3),a2->bondVectorTo(a1));}
	/**@brief Calculate the dihedral between any four atoms*/
	static inline double getDihedral(Atom* a1,Atom* a2,Atom* a3,Atom* a4) {return quickgeom::angleBetween(a1->bondVectorTo(a2).cross(a2->bondVectorTo(a3)),a2->bondVectorTo(a3).cross(a3->bondVectorTo(a4)));}
	/**@brief Calculate the improper dihedral between any four atoms. The central atom is first*/
	static inline double getImproper(Atom* a1,Atom* a2,Atom* a3,Atom* a4) {return quickgeom::angleBetween(a1->bondVectorTo(a2).cross(a1->bondVectorTo(a3)),a1->bondVectorTo(a3).cross(a1->bondVectorTo(a4)));}
}
#endif
