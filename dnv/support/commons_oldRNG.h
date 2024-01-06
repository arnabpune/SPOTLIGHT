#ifndef H_COMMONS
#define H_COMMONS 1
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <cmath>
#include<utility>
#include "../Eigen337/Eigen/Core"
#include "../Eigen337/Eigen/Dense"
#include "../Eigen337/Eigen/Geometry"
#define QUEUE_LIMIT 1000
using namespace std;

//Exceptions
/**This exception is thrown when any required file is not found*/
class FileNotFoundException : public exception {}; //Self-explanatory
/**This exception is thrown when any element, expected to be present in a list is not found in the list*/
class ElementNotFoundException : public exception {}; //Self-explanatory

static bool log_opened=false;
static const double PI=3.141592653589; /**<The constant PI with high accuracy*/
static double DIEL=138.0/10.0;/**<Dielectric constant used in simulation. This might need to be changed when a new force-field (other than CHARMM-27 is used)*/ //minus_one_by_beta=-1.0/beta;
static const double ANGVAR=(2.5*PI)/180.0; /**<The angle variance (gaussian error)*/
static const double DIHVAR=ANGVAR*0.75; /**<The dihedral variance (gaussian error)*/

//inline double sqr(double d) {return d*d;}
inline double sqr(const double& d) {return d*d;}

//Random number generation
/*namespace RNG
{
	static std::vector<mt19937> randomNumberGenerators;
}*/
/**@brief Get a random number in the uniform real interval (min,max)*/
double throwarandompoint(double min=0, double max=1){
	random_device rd; mt19937 eng(rd());
	uniform_real_distribution<> distr(min,max);
	return distr(eng);
}
/**@brief Get a random number in the normal distribution with given mean and standard deviation*/
double throwarandompoint_normal(double mean=0,double stdev=1.0)
{
	random_device rd; mt19937 eng(rd());
	normal_distribution<double> distr(mean,stdev);
  return distr(eng);
}
/**@brief Gives "true" or "false" as output with equal chance*/
inline bool toss() {return throwarandompoint(0,1)<0.5;}

//Template function(s)
template<class T,class U> const U& getCorresponding(const std::vector<std::pair<T,U>>& v,const T& t)
{
  for(auto& a : v)
  {
    if((get<0>(a))==t) return (get<1>(a));
  }
  cout << "Match not found\n";
  //return U();
}
/**@brief Select a random element from a list where each item has equal chance of being chosen*/
template<class T> T& randomSelect(std::vector<T>& v) {return v[(int)(throwarandompoint(0,v.size()))];}
/**@brief Select a random element from a list where each item has equal chance of being chosen. An alternative function for a "const" vector*/
template<class T> const T& randomSelect(const std::vector<T>& v) {return v[(int)(throwarandompoint(0,v.size()))];}
template<class T> std::ostream& operator<<(std::ostream& os,const std::vector<T>& v)
{
	os << "{";
	for(const T& t : v) os << t<<",";
	os << "}";
	return os;
}
//Some common methods
template<class T> inline void swap(std::vector<T>& v,int i,int j) {T temp=v[i]; v[i]=v[j]; v[j]=temp;}
/**@brief Sort a list using a sorting function that takes elements of the list as input and gives a comparable data-type (such as int or double) as output*/
template<class T,class U> void sort(std::vector<T>& v,U (*sortf)(const T& t))
{
	if(!v.size()) return;
	U min=sortf(v[0]);
	int is=0;
	for(int i=0;i<v.size();i++)
	{
		min=sortf(v[i]);
		is=i;
		for(int j=i+1;j<v.size();j++)
		{
			U temp=sortf(v[j]);
			if(temp<min) {min=temp; is=j;}
		}
		if(is!=i) swap(v,i,is);
	}
}
/**@brief Get the index of a particular element in a list*/
template<class T> int indexOf(const T& t,std::vector<T> ta)
{
	for(int i=0;i<ta.size();i++)
	{
		if(ta[i]==t) return i;
	}
	return -1;
}
template<class T> T* copyArray(const T* ptr,int n)
{
	T* cp=new T[n];
	for(int i=0;i<n;i++) cp[i]=ptr[i];
	return cp;
}
template<class T> std::vector<T> schuffle(std::vector<T> arr)
{
	std::vector<T> ret;
	int ind;
	while(arr.size())
	{
		ind=(int)throwarandompoint(0,ret.size());
		ret.push_back(arr[ind]);
		arr.erase(arr.begin()+ind);
	}
	return ret;
}
template<class T,class U> inline bool contains(std::pair<T,U> p,const T& o) {return ((get<0>(p))==o || (get<1>(p))==o);}
template<class T,class U> inline bool contains(std::pair<T,U> p,const U& o) {return ((get<0>(p))==o || (get<1>(p))==o);}
template<class T> inline bool contains(const std::vector<T>& v,const T& t)
{
	for(const T& el : v) {if(el==t) return true;}
	return false;
}
/**@brief Get the second element in a pair, given one*/
template<class T,class U> inline U getSecond(std::pair<T,U> p,const T& t)
{
	if((get<0>(p))==t) return (get<1>(p));
	else if((get<1>(p))==t) return (get<0>(p));
	else {cout <<"Warning: None in pair (getSecond: ~Line 220 in Molecule.cpp)";}
}
template<class T,class U> U findByKey(const std::vector<std::pair<T,U>>& v,const T& k)
{
	for(const auto& p : v) {if(get<0>(p)==k) return get<1>(p);}
	throw ElementNotFoundException();
}
inline static double dot_product(double xa,double ya,double za,double xb,double yb,double zb) {return xa*xb+ya*yb+za*zb;}
inline static double* cross_product(double xa,double ya,double za,double xb,double yb,double zb)
{
	double *cross=new double[3];
	cross[0]=ya*zb-za*yb;
	cross[1]=xb*za-xa*zb;
	cross[2]=xa*yb-ya*xb;
	return cross;
}
ofstream& fileExists(const std::string& logfname,char mode)
{
	char acc='N';
	if(std::ifstream(logfname) && !log_opened)
	{
		cout << "File ("<<logfname<<") exists. Append(A)/Overwrite(O)/Rename(R)/Skip(S)/Halt(H)?: ";
		cin >> acc;
		cout << acc<<"\n";
		if(acc=='H' || acc=='h') {cout << "User interrupt.\n"; exit(1);}
		else if(acc=='A' || acc=='a') {} //Do nothing
		else if(acc=='O' || acc=='o')  {auto tfs=ofstream(logfname,fstream::out); tfs.close();} //Overwritten
		else if(acc=='S' || acc=='s') {cout <<"Skipping logfile generation! Note: You will NOT have ANY DATA of the generated molecules except the final structure (unoptimized)\nType any letter to continue (or ^C to quit)..."; cin >>acc; }
		else if(acc=='R' || acc=='r')
		{
			cout << "Enter new name for logfile (no spaces or special characters other than _ and -): ";
			std::string fn;
			cin >> fn;
			return fileExists(fn,mode);
		}
		else {cout << "Input not understood. Type only the first character and try again.\n"; return fileExists(logfname,mode);}
	}
	ofstream& of=*(new ofstream()); if(acc!='S') {if(mode=='a') of.open(logfname,fstream::app); else of.open(logfname,fstream::out); log_opened=true;}
	if(of) return of;
	else {cout << "Log could not be opened\n"; return of;}
}

/**
	This namespace encompasses many of the string functions used in processing the data files
*/
namespace stringfx
{
	/**@brief Empty all remaining data from a stringstream into a string*/
	std::string drain(stringstream& ss)
	{
		std::string ret="",temp;
		while(true)
		{
			ss >> temp;
			ret+=" "+temp;
			if(ss.tellg()==-1) break;
		}
		return ret.substr(1,ret.size()-1);
	}
	std::string drain(istream& is)
	{
		std::string ret="",temp;
		while(true)
		{
			getline(is,temp);
			if(is.eof()) break;
			ret+=temp+"\n";
		}
		return ret;
	}
	inline int indexOf(char c,const std::string& s) {return s.find(c);}
	inline int indexOf(const char* st,const std::string& s) {return s.find(st);}
	inline int indexOf(const std::string& st,const std::string& s) {return s.find(st);}
	/**@brief split a string by a given delimiter and generate a std::vector of strings as output*/
	std::vector<std::string> split(const std::string& s,char delim=' ',bool autotrim=false)
	{
		std::vector<std::string> ret;
		std::string temp="";
		for(int i=0;i<s.length();i++)
		{
			if(s[i]==delim)
			{
				if(!temp.length() && autotrim) continue;
				ret.push_back(temp); temp=""; continue;
			}
			temp+=s[i];
		}
		ret.push_back(temp);
		return ret;
	}
	inline std::string& replace(std::string& s,std::string val,std::string ns)
	{
		int l=s.find(val);
		if(l==std::string::npos) return s;
		return s.replace(l,val.size(),ns);
	}
	inline std::string trim(const std::string& aString) {if(aString.find_first_not_of(' ')==std::string::npos) return ""; return aString.substr(aString.find_first_not_of(' '), (aString.find_last_not_of(' ') - aString.find_first_not_of(' ')) + 1);}

	bool matches(const std::string& seq,const std::string& str)
	{
		if(seq==str) return true;
		bool wild=false,pos=false;
		std::string cmp=seq;
		if(seq[seq.size()-1]=='*') {wild=true; pos=0; cmp=seq.substr(0,seq.size()-1);}
		if(seq[0]=='*') {wild=true; pos=1; cmp=seq.substr(1,seq.size()-1);}
		if(!wild) return false;
		if(!pos) return (indexOf(cmp,str)==0);
		else return (indexOf(cmp,str)==str.size()-cmp.size());
	}
}
namespace quickmath
{
	std::vector<double> getGradient(double (*func)(std::vector<double>),std::vector<double> point,double D=1e-6)
	{
		if(point.size()<1) return std::vector<double>();
		std::vector<double> ret,newpoint;
		double df;
		for(int i=0;i<point.size();i++)
		{
			newpoint=std::vector<double>();
			for(int j=0;j<i;j++) newpoint.push_back(point[j]);
			newpoint.push_back(point[i]+D);
			for(int j=i+1;j<point.size();j++) newpoint.push_back(point[j]);
			df=func(newpoint)-func(point);
			ret.push_back(df/D);
		}
		return ret;
	}
	std::vector<double> gradientDescentFit(double (*func)(std::vector<double>),int size,double err=1e-5,std::vector<double> initguess=std::vector<double>()) //Func is the function to minimize
	{
		if(size<1) return std::vector<double>();
		if(initguess.size()!=size) {initguess=std::vector<double>(); for(int i=0;i<size;i++) initguess.push_back(1);}
		std::vector<double> guess=initguess;
		int K=0;
		double D=err*10;
		double ov=1e5,tv;
		bool acc=true;
		cout << "D="<<D<<"\n";
		while(abs(func(guess)-ov)>err)
		{
			//if(K++>10000) {cout <<"Turn limit reached!"; break;}
			cout << guess << "\n";
			ov=func(guess);
			std::vector<double> grad=getGradient(func,guess,err);
			//acc=false;
			for(int i=0;i<size;i++)
			{
				tv=grad[i]*D;
				cout << tv << ",";
				//if(abs(tv)<err) continue;
				guess[i]-=tv;
				//acc=true;
			}
			cout << "\n";
			//if(!acc) break;
			if(guess[0]!=guess[0]) break;
		}
		cout << "Residual: "<<func(guess)<<"\n";
		return guess;
	}
}
/**
	This namespace contains all the methods used for geometric evaluations in DeNovo. This includes position trialling, angle calculation, and random radial vector generation.<br/>Note: All angle related methods take angle input in radians and output is given also in radians
	This namespace makes extensive use of the Eigen libraries. See Eigen: http://eigen.tuxfamily.org/ for more details. (git repo: https://gitlab.com/libeigen/eigen.git)
*/
namespace quickgeom
{
	/**@brief Get literally <b>any</b> random unit vector*/
	static inline Eigen::Vector3d getRandomUnitVector()
	{
	  double theta,phi;
	  theta=throwarandompoint(0,PI); phi=2*throwarandompoint(0,PI);
	  return Eigen::Vector3d(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
	}
	/**@brief Get a random vector pointing to a point at a given distance from a fixed center
		 @details This method generates a vector pointing to a point dispersed randomly around a fixed center at a given distance from it
		 @param[in] origin: The center about which the points are to be dispersed
		 @param[in] radius: The distance from the center for all the points
	*/
	static inline Eigen::Vector3d getRandomVectorInShell(const Eigen::Vector3d& origin,const double& radius) {return radius*getRandomUnitVector()+origin;}
	/**@brief A wrapper to generate multiple vectors around a given centre (useful while trialling seed positions or the position of the second atom). See getRandomVectorInShell(const Eigen::Vector3d&,const double&)*/
	static std::vector<Eigen::Vector3d> getRandomVectorsInShell(const Eigen::Vector3d& origin,const double& radius,int sampleSize)
	{
	  std::vector<Eigen::Vector3d> ret;
	  for(int i=0;i<sampleSize;i++)
	    ret.push_back(getRandomVectorInShell(origin,radius));
	  return ret;
	}
	/**@brief Get the angle subtended by two vectors*/
	static inline double angleBetween(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2) {return acos(v1.dot(v2)/(v1.norm()*v2.norm()));}
	/**@brief Get the normal unit vector to the plane created by two vectors*/
	static inline Eigen::Vector3d planeOf(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2) {return (v1.cross(v2))/(v1.norm()*v2.norm());}
	/**@brief Given three bond vectors, calculate the dihedral angle they make*/
	static inline double dihedralBetween(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2, const Eigen::Vector3d& v3) {return angleBetween(v1.cross(v2),v2.cross(v3));}

	/**Get a rotation matrix for rotation about a particular axis, by a given angle*/
	static inline auto getRotationMatrix(const Eigen::Vector3d& axis,double angle) {return Eigen::AngleAxisd(angle,axis);} //angle in radians
	static std::vector<Eigen::Vector3d> getConePositions(const Eigen::Vector3d& fixed,double length,double angle,int num=10) //angle in radians
	{
		Eigen::Vector3d axis=Eigen::Vector3d::UnitZ().cross(fixed/fixed.norm());
		auto rM=getRotationMatrix(axis/axis.norm(),acos(fixed(2)/fixed.norm()));
		//rM= pow(Eigen::Matrix<double,3,3>(rM).determinant(),1/3.0)*rM;
		std::vector<Eigen::Vector3d> poses;
		double phi;
		double theta;
		for(int i=0;i<num;i++)
		{
			phi=throwarandompoint(0,2*PI);
			theta=throwarandompoint_normal(angle,ANGVAR);
			poses.push_back(length*(rM*Eigen::Vector3d(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))));
		}
		return poses;
	}
	/*static std::vector<Eigen::Vector3d> getArcPositions(const Eigen::Vector3d& fixed,const Eigen::Vector3d& def,double length,double angle,std::vector<double> dihedral,int num=10)
	{
		Eigen::Vector3d plane=fixed.cross(def); plane/=plane.norm();
		Eigen::Vector3d axisd=Eigen::Vector3d::UnitY().cross(plane/plane.norm());
		auto nrM=getRotationMatrix(axisd/axisd.norm(),acos(plane(1)/plane.norm()));
		Eigen::Vector3d nfixed=fixed;
		Eigen::Vector3d axis=Eigen::Vector3d::UnitZ().cross(nfixed/nfixed.norm());
		//if(angleBetween(plane,Eigen::Vector3d::UnitY())>angleBetween(plane,-Eigen::Vector3d::UnitY())) plane=-plane;
		auto rM=getRotationMatrix(axis/axis.norm(),acos(nfixed(2)/nfixed.norm()));
		cout <<"Plane angle: "<< angleBetween(Eigen::Vector3d::UnitY(),plane) << "\n";
		//auto nrM=getRotationMatrix(fixed/fixed.norm(),acos(plane(1)/plane.norm()));
		rM=nrM;
		//cout <<"Rotation determinant: "<< rM.determinant() << "\n";
		//double na=angleBetween((rM*Eigen::Vector3d(sin(angle),0,cos(angle))).cross(fixed),plane);
		//for(double& d : dihedral) d+=na;
		std::vector<Eigen::Vector3d> poses;
		double phi;
		double theta;
		double dih;
		for(int i=0;i<num;i++)
		{
			theta=throwarandompoint_normal(angle,ANGVAR);
			if(dihedral.size()==1) phi=throwarandompoint_normal(dihedral[0],DIHVAR);
			else phi=throwarandompoint_normal(randomSelect(dihedral),DIHVAR);
			cout << theta <<","<<sin(theta)<<"\t"<<phi<<","<<sin(phi)<<"\n";
			poses.push_back(length*(rM*Eigen::Vector3d(sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta))));
		}
		return poses;
	}*/
	/**@brief The function used for generating trial positions with a given length, angle and set of allowed dihedrals
		 @details This function uses geometric calculations to choose trial positions when placing a new atom. It also adds in the gaussian error. See the algorithm PDF for details on the calculation.
		 @param[in] fixed: The primary axis vector (The vector representing the bond whose one end will be connected  to the atom being trialled)
		 @param[in] def: The defining vector (or fixed vector). It is the bond-vector with respect to which the dihedrals are defined.
		 @param[in] length: The bond length. The vectors returned will have this magnitude, ensuring that the new atom is placed at this distance from the source.
		 @param[in] angle: The bond angle. This angle will be maintained from the "fixed" vector. See algorithm PDF for details on the exact calculation.
		 @param[in] dihedral: The set of all allowed dihedrals in a vector object
		 @param[in] num: Number of trial positions to generate
		 @return A list of Vectors (Eigen::Vector3d) which has the suggested positions. All vectors are centered at origin and have to be translated to the source atom.
	*/
	static std::vector<Eigen::Vector3d> getArcPositions(const Eigen::Vector3d& fixed,const Eigen::Vector3d& def,double length,double angle,const std::vector<double>& dihedral,int num=10)
	{
		double beta=fixed.dot(def)/(fixed.norm()*def.norm());
		//Eigen::Vector3d plane=fixed.cross(def); plane/=plane.norm();
		double theta,phi;
		std::vector<Eigen::Vector3d> poses;
		Eigen::Vector3d prpos,per=def-(def.dot(fixed)/fixed.norm())*(fixed/fixed.norm()); per/=per.norm(); //per=def/def.norm()-(def.dot(fixed/fixed.norm())*fixed/fixed.norm()); per/=per.norm();
		for(int i=0;i<num;i++)
		{
			theta=throwarandompoint_normal(angle,ANGVAR);
			if(dihedral.size()==1) phi=throwarandompoint_normal(dihedral[0],DIHVAR);
			else phi=throwarandompoint_normal(randomSelect(dihedral),DIHVAR);
			//double x=((toss())?1:-1)*acos((cos(theta)*beta - sqrt(beta*beta + sin(theta)*sin(theta)))/(1+beta*beta));
			//prpos=(sin(x)*(fixed/fixed.norm())+cos(x)*(def/def.norm()));
			prpos=cos(theta)*(fixed/fixed.norm())+sin(theta)*per;
			//cout << prpos.norm()<<",";
			//if(prpos.cross(fixed).dot(plane)>0) phi+=PI;
			auto nrM=getRotationMatrix(fixed/fixed.norm(),(toss()?-1:1)*phi);
			poses.push_back(nrM*(length*prpos));
		}
		//cout << "\n";
		return poses;
	}
	/**@brief A wrapper for calling getArcPositions() with a single allowed dihedral instead of a set*/
	inline static std::vector<Eigen::Vector3d> getArcPositions(const Eigen::Vector3d& fixed,const Eigen::Vector3d& def,double length,double angle,double dihedral,int num=10)
	{
		std::vector<double> dihs(1,dihedral); //dihs.push_back(dihedral);
		return getArcPositions(fixed,def,length,angle,dihs,num);
	}

	//Newer methods added in after THETA2
	class Cuboid;
	class Sphere;
	/**\class Container commons.h "support/commons.h"
		 This is the (abstract) Container class<br/>
		 It allows the user to design many shapes of containers as long as distances to various points can be computed as required by the virtual functions (see functions of this class)<br/>
		 The initial intent of this class is to provide a container to restrain molecule growth (as an alternative to using restraints to hotspot residues)<br/>
		 See Also: Cuboid, Sphere
	*/
	class Container
	{
	protected:
		Eigen::Vector3d centre; /**<The geometric centre for this container*/
		/**@brief Default empty constructor. Takes the centre to be at (0,0,0)*/
		Container() : Container(0,0,0) {}
		/**@brief Construct an abstract container with centre at a point given by the inpue vector*/
		Container(const Eigen::Vector3d& v) {centre=v;}
		/**@brief Construct an abstract container with centre at (x,y,z)*/
		Container(double x,double y, double z) : Container(Eigen::Vector3d(x,y,z)) {}

	public:
		/**@brief Get the geometric centre of this container (an alternative for a "const" Container object)*/
		inline const Eigen::Vector3d& getCentre() const {return centre;}
		/**@brief Get the geometric centre of this container*/
		inline Eigen::Vector3d getCentre() {return centre;}

		/**@brief Checks if the given point is present inside this region. The boundary points are not included*/
		virtual inline bool contains(double x, double y, double z) const =0;
		/**@brief Checks if the given point is present inside this region. The boundary points are not included*/
		virtual inline bool contains(const Eigen::Vector3d& r) const {return contains(r(0),r(1),r(2));}
		/**@brief Get the distance of this point to the centre of this region*/
		inline double getDistanceFrom(double x, double y, double z) const {return sqrt(sqr(centre(0)-x)+sqr(centre(1)-y)+sqr(centre(2)-z));}
		/**@brief Get the distance of this point to the centre of this region*/
		inline double getDistanceFrom(const Eigen::Vector3d& r) const {return getDistanceFrom(r(0),r(1),r(2));}
		/**@brief Get the distance from the given point to the closest point on the surface of the region. It is zero if contained inside the region*/
		virtual inline double getMinimumDistanceFrom(const Eigen::Vector3d& v) const {return getMinimumDistanceFrom(v(0),v(1),v(2));}
		/**@brief Same as getMinimumDistanceFrom(const Eigen::Vector3d&) but explicitly specifies the coordinates*/
		virtual inline double getMinimumDistanceFrom(double x, double y, double z) const =0;
		/**@brief Get the volume of this container*/
		virtual inline double getVolume() const = 0;

		/**@brief Get the smallest possible cuboid enclosing this region*/
		virtual inline const Cuboid getContainingBox() const =0;
		/**@brief Get the smallest possible sphere enclosing this region*/
		virtual inline const Sphere getContainingSphere() const =0;
	};

	class Cuboid : public Container
	{
		Eigen::Vector3d dim;
	public:
		Cuboid() : Cuboid(0,0,0,1,1,1) {}
		/**@brief Construct a cubic container using two vectors
			 @details The vectors specify the box parameters
			 @param[in] dimen: The dimensions (size) as (x,y,z)
			 @param[in] v: The positional vector parameter (see the 'center' parameter)
			 @param[in] center: If set to true, 'v' represents the position of center of the cuboid, otherwise, it denotes the origin (bottom-left-front corner) (Default - v is centre)
		*/
		Cuboid(const Eigen::Vector3d& v,const Eigen::Vector3d& dimen,bool center=true) : Container()
		{
			dim=dimen;
			if(center) centre=v;
			else centre=Eigen::Vector3d(v(0)+dimen(0)/2.0,v(1)+dimen(1)/2.0,v(2)+dimen(2)/2.0);
		}
		/**@brief Same as Cuboid(const Eigen::Vector3d&,const Eigen::Vector3d&,bool) but explicitly uses doubles instead of Eigen::Vector3d objects*/
		Cuboid(double cx,double cy,double cz,double sx,double sy,double sz,bool center=true) : Cuboid(Eigen::Vector3d(cx,cy,cz),Eigen::Vector3d(sx,sy,sz),center) {}
		/**@brief Same as Cuboid(const Eigen::Vector3d&,const Eigen::Vector3d&,bool) but uses mixed data-types instead of just Eigen::Vector3d objects*/
		Cuboid(const Eigen::Vector3d& v,double sx,double sy,double sz,bool center=true) : Cuboid(v,Eigen::Vector3d(sx,sy,sz),center) {}
		/**@brief Same as Cuboid(const Eigen::Vector3d&,const Eigen::Vector3d&,bool) but uses mixed data-types instead of just Eigen::Vector3d objects*/
		Cuboid(double cx,double cy,double cz,const Eigen::Vector3d& dimen,bool center=true) : Cuboid(Eigen::Vector3d(cx,cy,cz),dimen,center) {}

		inline bool contains(double x,double y,double z) const override {return ((x<centre(0)+dim(0)/2 && x>centre(0)-dim(0)/2) && (y<centre(1)+dim(1)/2 && y>centre(1)-dim(1)/2) && (z<centre(2)+dim(2)/2 && z>centre(2)-dim(2)/2));}
		inline double getMinimumDistanceFrom(double x,double y,double z) const override
		{
			double xdist=(centre(0)<x)?max(0.0,x-(centre(0)+dim(0)/2.0)):max(0.0,(centre(0)-dim(0)/2.0)-x);
			double ydist=(centre(1)<y)?max(0.0,y-(centre(1)+dim(1)/2.0)):max(0.0,(centre(1)-dim(1)/2.0)-y);
			double zdist=(centre(2)<z)?max(0.0,z-(centre(2)+dim(2)/2.0)):max(0.0,(centre(2)-dim(2)/2.0)-z);
			return sqrt(sqr(xdist)+sqr(ydist)+sqr(zdist));
		}

		inline const Cuboid getContainingBox() const override {return *this;}
		inline const Sphere getContainingSphere() const override;

		inline double getVolume() const override {return dim(0)*dim(1)*dim(2);}
		inline double getMinX() const {return centre(0)-dim(0)/2;}
		inline double getMaxX() const {return centre(0)+dim(0)/2;}
		inline double getMinY() const {return centre(1)-dim(1)/2;}
		inline double getMaxY() const {return centre(1)+dim(1)/2;}
		inline double getMinZ() const {return centre(2)-dim(2)/2;}
		inline double getMaxZ() const {return centre(2)+dim(2)/2;}
	};
	class Sphere : public Container
	{
		double rad; /**<Radius of sphere*/
	public:
		/**@brief Default constructor that constructs a Sphere of radius 1nm around (0,0,0)*/
		Sphere() : Sphere(0,0,0,1) {}
		/**@brief Construct a sphere around a given point with a given radius*/
		Sphere(const Eigen::Vector3d& c,double r)
		{
			centre=c;
			rad=r;
		}
		/**@brief Same as Sphere(const Eigen::Vector3d&,double) except that the centre is explicitly mentioned*/
		Sphere(double cx,double cy,double cz,double r) : Sphere(Eigen::Vector3d(cx,cy,cz),r) {}

		inline bool contains(double x,double y,double z) const override {return (this->getDistanceFrom(x,y,z)<rad);}
		inline double getMinimumDistanceFrom(double x,double y,double z) const override {return max(this->getDistanceFrom(x,y,z)-rad,0.0);}

		inline const Sphere getContainingSphere() const override {return *this;}
		inline const Cuboid getContainingBox() const override;

		inline double getVolume() const override {return (4/3.0)*PI*rad*rad*rad;}
	};

	inline const Sphere Cuboid::getContainingSphere() const {return Sphere(centre,dim.norm()/2.0);} //To be set once Sphere class is defined
	inline const Cuboid Sphere::getContainingBox() const {return Cuboid(centre,rad*2,rad*2,rad*2,true);}

	/** \class ContainerUnion commons.h "support/commons.h"
		 The ContainerUnion class allows to take the 'union' of multiple containers as a region. Any point contained in ANY ONE of the containers will be in the union.<br/>
		 The containing box (See Container::getContainingBox()) can sphere (See Container::getContainingSphere()) will contain the entire region.<br/>
	*/
	class ContainerUnion : public Container
	{
		std::vector<const Container*> containers;
	public:
		/**@brief Construct a unified region which is the union of all the containers (supplied as pointers) provided in the list*/
		ContainerUnion(const std::vector<const Container*>& conts)
		{
			assert(conts.size());
			containers=conts;
			Eigen::Vector3d nc;
			double sumvols=0;
			for(const Container* cnt : conts)
			{
				double vol=cnt->getVolume();
				nc+=cnt->getCentre()*vol;
				sumvols+=vol;
			}
			centre=nc/sumvols;
		}
		/**@brief Construct an union with a single onject*/
		ContainerUnion(const Container* single) : ContainerUnion(std::vector<const Container*>(1,single)) {}

		inline bool contains(double x,double y,double z) const override
		{
			for(const Container* cnt : containers) {if(cnt->contains(x,y,z)) return true;}
			return false;
		}
		inline double getMinimumDistanceFrom(double x,double y,double z) const override
		{
			double mindist=containers[0]->getMinimumDistanceFrom(x,y,z);
			double tempd;
			int i=1;
			while(mindist!=0 && i<containers.size()) {if((tempd=containers[i]->getMinimumDistanceFrom(x,y,z))<mindist) mindist=tempd; i++;}
			return mindist;
		}
		inline const Cuboid getContainingBox() const override
		{
			const Cuboid& c1=containers[0]->getContainingBox();
			double minx=c1.getMinX(),miny=c1.getMinY(),minz=c1.getMinZ(),maxx=c1.getMaxX(),maxy=c1.getMaxY(),maxz=c1.getMaxZ();
			for(int i=1;i<containers.size();i++)
			{
				const Cuboid& tempcbd=containers[i]->getContainingBox();
				if(minx>tempcbd.getMinX()) minx=tempcbd.getMinX();
				if(miny>tempcbd.getMinY()) miny=tempcbd.getMinY();
				if(minz>tempcbd.getMinZ()) minz=tempcbd.getMinZ();
				if(maxx<tempcbd.getMaxX()) maxx=tempcbd.getMaxX();
				if(maxy<tempcbd.getMaxY()) maxy=tempcbd.getMaxY();
				if(maxz<tempcbd.getMaxZ()) maxz=tempcbd.getMaxZ();
			}
			return Cuboid(minx,miny,minz,maxx,maxy,maxz,false);
		}
		inline const Sphere getContainingSphere() const override {return this->getContainingBox().getContainingSphere();}
		inline double getVolume() const override
		{
			double v;
			for(const Container* cnt : containers) v+=cnt->getVolume();
			return v;
		}
	};
}
/** \class Stack commons.h "support/commons.h"
	 This class is DeNovo's implementation of the standard stack with a constant maximum limit of allowed elements.<br/>
	 It will not be explained in detail
*/
template<class T> struct Stack
{
	T* head;
	T* tail;
	int size;
	Stack(int size)
	{
		tail=new T[size];
		head=tail;
	}
	~Stack() {delete[] tail;}

	/*@brief Get the maximum number of elements this stack can hold*/
	inline int getMaxSize() const {return size;}
	/**@brief Get the number of elements in this stack*/
	inline int getSize() const {return (head-tail);}
	/**@brief Get the next element (and remove it from the stack)*/
	inline T next() {return *(--head);}
	/**@brief Add an element to the stack*/
	inline T* push(const T& obj) {*head=obj; return head++;}
	/**@brief Is the stack empty?*/
	inline bool isEmpty() const {return (tail==head);}
	/**@brief Empty the stack into a vector and return it.*/
	std::vector<T> drain()
	{
		std::vector<T> ret;
		while(!isEmpty()) ret.push_back(next());
		return ret;
	}
};
/** \class Queue commons.h "support/commons.h"
	 This class is DeNovo's implementation of the standard cyclic queue with a constant maximum limit of allowed elements.<br/>
	 It will not be explained in detail
*/
template<class T> struct Queue
{
	T* start;
	T* head;
	T* tail;
	int size;
	Queue(int s)
	{
		size=s;
		start=new T[size];
		tail=start;
		head=tail;
	}
	~Queue() {delete[] start;}

	/*@brief Get the maximum number of active elements this queue can hold.*/
	inline int getMaxSize() const volatile {return size;}
	/**@brief Get the number of elements actively in this queue*/
	inline int getSize() const volatile {return (head>=tail)?(head-tail):(size-(tail-head));}
	/**@brief Get the next element (and remove it from queue)*/
	inline T next() volatile {if(tail-start>=size) tail=start; return *(tail++);}
	/**@brief Add an element to the queue*/
	inline T* push(const T& obj) volatile {if(head-start>=size) head=start; *head=obj; return head++;}
	/**@brief Is the queue empty?*/
	inline bool isEmpty() const volatile {return (tail==head);}
	/**@brief Empty the queue into a vector and return it.*/
	std::vector<T> drain() volatile
	{
		std::vector<T> ret;
		while(!isEmpty()) ret.push_back(next());
		return ret;
	}
};
#endif
