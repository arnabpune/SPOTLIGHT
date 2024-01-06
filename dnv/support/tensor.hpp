/*
We will first make a tensor object and then segment it to double and int types
*/

#include "./boxgrid_based_operations.hpp"
#pragma once
using namespace std;

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#######################Define a 3-dimensional tensor########################|
|#########################for storing DOUBLES only###########################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/

class tensor3D_double{
	private:
		vector<Eigen::MatrixXd> layers;
	public:
		int wid,hei,dep;
		//only a default constructor is needed:
		tensor3D_double() {wid=hei=dep=0;}
		tensor3D_double(int,int,int);
		tensor3D_double(const tensor3D_double& o);
		~tensor3D_double();
		//Define some filling functions
		void fill_with_zero();
		void fill_with_identity();
		void fill_with_one();
		void fill_with_number(double);
		void fill_with_random();
		//Modification of elements
		void modify(int,int,int,double);
		//Accessing elements
		double &operator()(int i,int j,int k){
			Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,dep,hei);
			int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
			Eigen::MatrixXd m;
			m = layers[k_c];
			return m(i_c,j_c);
		}
		const double &operator()(int i,int j,int k) const{
			Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,dep,hei);
			int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
			Eigen::MatrixXd m;
			m = layers[k_c];
			return m(i_c,j_c);
		}
		double get(int,int,int) const;
		Eigen::MatrixXd get_layer(int) const;
};

/*
We are going to define the constructors, copy-constructor and destructor
*/

tensor3D_double::tensor3D_double(int w,int h,int d){
	layers.reserve(d);
	int i = 0;
	while(i < d){
		Eigen::MatrixXd m(w,h);
		layers.push_back(m);
		i = i + 1;
	}
	wid=w; hei=h; dep=d;
}

tensor3D_double::tensor3D_double(const tensor3D_double& o){
	vector<Eigen::MatrixXd>::iterator it; it = layers.begin();
	int i = 0;
	while(it != layers.end()){
		*it = o.get_layer(i);
		i = i + 1; ++it;
	}
	wid=o.wid; hei=o.hei; dep=o.dep;
}

tensor3D_double::~tensor3D_double(){}

/*
We will define the filler functions so that we can define constant matrices or random matrices
*/
inline void tensor3D_double::fill_with_zero() {fill_with_number(0);}

void tensor3D_double::fill_with_identity(){
	vector<Eigen::MatrixXd>::iterator it; it = layers.begin();
	while(it != layers.end()){
		*it = Eigen::MatrixXd::Identity(wid,hei);
		++it;
	}
}

inline void tensor3D_double::fill_with_one(){ this->fill_with_number(1);}

void tensor3D_double::fill_with_number(double val){
	vector<Eigen::MatrixXd>::iterator it; it = layers.begin();
	while(it != layers.end()){
		*it = Eigen::MatrixXd::Constant(wid,hei,val);
		++it;
	}
}

void tensor3D_double::fill_with_random(){
	vector<Eigen::MatrixXd>::iterator it; it = layers.begin();
	while(it != layers.end()){
		*it = Eigen::MatrixXd::Random(wid,hei);
		++it;
	}
}

/*
We will define a modifier function to modify a specific element of the tensor
*/

void tensor3D_double::modify(int i,int j,int k,double d){
	Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,hei,dep);
	int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
	vector<Eigen::MatrixXd>::iterator it; it = layers.begin();
	int kk = 0;
	vector<Eigen::MatrixXd> layer_temp; layer_temp.reserve(hei);
	while(it != layers.end()){
		if( (kk = k_c) ){
			//this is the layer we are looking for
			Eigen::MatrixXd m(wid,hei); m = *it; //assign the matrix
			m(i_c,j_c) = d;//modify the element
			layer_temp.push_back(m);
			kk++; ++it;
		}
		else{
			Eigen::MatrixXd m(wid,hei); m = *it;
			layer_temp.push_back(m);
			kk++; ++it;
		}
	}
	layers = layer_temp;
}

/*
We define an operator overload to access elements in the object like an array
using general mathematical indexing starting from 0 ofcourse, this is then repeated in the
get() function.
*/

double tensor3D_double::get(int i,int j,int k) const{
	Eigen::MatrixXd m(wid,hei);
	m = layers[k];
	double d = m(i,j);
	return d;
}

/*
We define a function to return a single layer or slice along the z-axis
*/
Eigen::MatrixXd tensor3D_double::get_layer(int k) const{
	Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(0,0,k,wid,hei,dep);
	int k_c = corrected_vector(2);
	Eigen::MatrixXd m(wid,hei);
	m = layers[k_c];
	return m;
}

// Overriding the << operator to print the tensor
std::ostream& operator<<(ostream& os, const tensor3D_double& t)
{
	for(int i = 0; i < t.wid; i++){
		for(int j = 0; j < t.hei; j++){
			for(int k = 0; k < t.dep; k++){
				os << i << " " << j <<" " << k<<" "<<t.get(i,j,k) << "\n";
			}
		}
	}
	return os;
}



/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#######################Define a 3-dimensional tensor########################|
|#########################for storing INTEGERS only##########################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/

//template<int Nx,int Ny,int Nz>
class tensor3D_int{
	private:
		vector<Eigen::MatrixXi> layers;
	public:
		int wid,hei,dep;
	public:
		//only a default constructor is needed:
		tensor3D_int() {wid=hei=dep=0;}
		tensor3D_int(int xw,int yw,int zw);
		tensor3D_int(const tensor3D_int& o);
		~tensor3D_int();
		//Define some filling functions
		void fill_with_zero();
		void fill_with_identity();
		void fill_with_one();
		void fill_with_number(int);
		void fill_with_random();
		//Modification of elements
		void modify(int,int,int,int);
		//Accessing elements
		int &operator()(int i,int j,int k){
			Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,hei,dep);
			int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
			return layers[k_c](i_c,j_c);
		}
		const int &operator()(int i,int j,int k) const{
			Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,hei,dep);
			int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
			return layers[k_c](i_c,j_c);
		}
		int get(int,int,int) const;
		Eigen::MatrixXi get_layer(int) const;
		void printTable() const;
		//void render();
};

/*
We are going to define the constructors, copy-constructor and destructor
*/
//template<int Nx,int Ny,int Nz>
tensor3D_int::tensor3D_int(int x,int y,int z){
	layers.reserve(z);
	int i = 0;
	while(i < z){
		Eigen::MatrixXi m(x,y);
		layers.push_back(m);
		i = i + 1;
	}
	wid=x; hei=y; dep=z;
}

tensor3D_int::tensor3D_int(const tensor3D_int& o){
	vector<Eigen::MatrixXi>::iterator it; it = layers.begin();
	int i = 0;
	while(it != layers.end()){
		*it = o.get_layer(i);
		i = i + 1; ++it;
	}
	wid=o.wid; hei=o.hei; dep=o.dep;
}

tensor3D_int::~tensor3D_int(){}

/*
We will define the filler functions so that we can define constant matrices or random matrices
*/

inline void tensor3D_int::fill_with_zero(){fill_with_number(0);}

void tensor3D_int::fill_with_identity(){
	vector<Eigen::MatrixXi>::iterator it; it = layers.begin();
	while(it != layers.end()){
		*it = Eigen::MatrixXi::Identity(wid,hei);
		++it;
	}
}

inline void tensor3D_int::fill_with_one(){ fill_with_number(1); }

void tensor3D_int::fill_with_number(int val){
	for(int i=0;i<dep;i++)
		layers[i]=Eigen::MatrixXi::Constant(wid,hei,val);
	cout << "Grid has "<<dep<<" layers, each of dimension: "<<wid<<"x"<<hei<<"\n";
}

void tensor3D_int::fill_with_random(){
	vector<Eigen::MatrixXi>::iterator it; it = layers.begin();
	while(it != layers.end()){
		*it = Eigen::MatrixXi::Random(wid,hei);
		++it;
	}
}

/*
We will define a modifier function to modify a specific element of the tensor
*/

void tensor3D_int::modify(int i,int j,int k,int d){
	Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(i,j,k,wid,hei,dep);
	int i_c = corrected_vector(0); int j_c = corrected_vector(1); int k_c = corrected_vector(2);
	layers[k_c](i_c,j_c)=d;
	/*vector<Eigen::MatrixXi>::iterator it; it = layers.begin();
	int kk = 0;
	vector<Eigen::MatrixXi> layer_temp; layer_temp.reserve(dep);
	while(it != layers.end()){
		if( (kk = k_c) ){
			//this is the layer we are looking for
			Eigen::MatrixXi m(wid,hei); m = *it; //assign the matrix
			m(i_c,j_c) = d;//modify the element
			layer_temp.push_back(m);
			kk = kk + 1; ++it;
		}
		else{
			Eigen::MatrixXi m(wid,hei); m = *it;
			layer_temp.push_back(m);
			kk = kk + 1; ++it;
		}
	}
	layers = layer_temp;*/
}

/*
We define an operator overload to access elements in the object like an array
using general mathematical indexing starting from 0 ofcourse, this is then repeated in the
get() function.
*/

inline int tensor3D_int::get(int i,int j,int k) const { return this->operator()(i,j,k);}

/*
We define a function to return a single layer or slice along the z-axis
*/
Eigen::MatrixXi tensor3D_int::get_layer(int k) const {
	Eigen::Vector3i corrected_vector; corrected_vector = pbc_index_correction(0,0,k,wid,hei,dep);
	int k_c = corrected_vector(2);
	Eigen::MatrixXi m(wid,hei);
	m = layers[k_c];
	return m;
}
void tensor3D_int::printTable() const
{
	for(int k=0;k<dep;k++)
	{
		const Eigen::MatrixXi& m=layers[k];
		for(int i=0;i<wid;i++)
		{
			for(int j=0;j<hei;j++)
				cout << i <<" "<< j <<" "<< k <<" "<< m(i,j) << "\n";
		}
	}
}
