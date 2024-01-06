#include "../assign_position_to_box.hpp"

using namespace std;

int main(){
	//Units are in nm
	int maxr = 5.0;
	double rcut = 1.0;
	//create random position vectors of random sizes
	int Npos = 50;
	while(Npos < 1000){
		cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
		cout<<"|############################################################################|"<<endl;
		cout<<"Currently performing with N = "<<Npos<<" elements"<<endl;
		int n = 0;
		//first generate a vector of double positions
		vector<Eigen::Vector3d> allpos;
		while( n < Npos ){
			Eigen::Vector3d randpos; randpos = genrandomvec3d(0.0,5.0,0.0,5.0,0.0,5.0);
			allpos.push_back(randpos);
			n = n + 1;
		}
		//figure out the position vector size
		cout<<"Current Position Vector size = "<<allpos.size()<<endl;
		//create the box and print it's associated variables
		cout<<"Creating box now"<<endl;
		pair<Eigen::Vector4d,Eigen::Vector3i> p; p = define_bounding_box(allpos, rcut);
		Eigen::Vector4d r4d; r4d = get<0>(p);
		Eigen::Vector3i abc; abc = get<1>(p);
		Eigen::Vector3d r_o(r4d(0),r4d(1),r4d(2)); double d = r4d(3);
		//const int a,b,c;
		int a = abc(0); int b = abc(1); int c = abc(2);
		cout<<"Center is at = "<<r_o(0)<<","<<r_o(1)<<","<<r_o(2)<<endl;
		cout<<"d-spacing is = "<<d<<" nm"<<endl;
		cout<<"Box vectors are: a = "<<a<<"nm, b = "<<b<<"nm, c = "<<c<<"nm"<<endl;
		//create the subdivisions Nx,Ny,Nz with grid bins and print it as well
		const int Nx = 50; const int Ny = 50; const int Nz = 50;
		cout<<"Guessed grid-bins are determined as: Nx = "<<a/d<<", Ny = "<<b/d<<", Nz = "<<c/d<<endl;
		cout<<"User supplied grid-bins are: Nx = "<<Nx<<", Ny = "<<Ny<<", Nz = "<<Nz<<endl;
		//because the user has supplied some parameters thus we need to perform all operations with this in mind
		d = min( (((double)c)/((double)Nz)),(min( (((double)a)/((double)Nx)), (((double)b)/((double)Ny)) )) );
		//create the tensors and fill them with a preset zero and print that out as well.
		tensor3D_double<Nx,Ny,Nz> newdoublegrid; tensor3D_int<Nx,Ny,Nz> newintgrid;
		newdoublegrid.fill_with_zero(); newintgrid.fill_with_zero();
		cout<<"Double and Integer 2nd order tensor grids created!"<<endl;
		//Try to modify a sample element in the grids
		cout<<"OLD element M_double("<<Nx-2<<","<<Ny-2<<","<<Nz-2<<") = "<<newdoublegrid.get(Nx-2,Ny-2,Nz-2)<<endl;
		newdoublegrid.modify(Nx-2,Ny-2,Nz-2,128.0);
		cout<<"NEW element M_double("<<Nx-2<<","<<Ny-2<<","<<Nz-2<<") = "<<newdoublegrid.get(Nx-2,Ny-2,Nz-2)<<endl;
		cout<<"OLD element M_int("<<Nx-2<<","<<Ny-2<<","<<Nz-2<<") = "<<newintgrid.get(Nx-2,Ny-2,Nz-2)<<endl;
		newintgrid.modify(Nx-2,Ny-2,Nz-2,128);
		cout<<"NEW element M_int("<<Nx-2<<","<<Ny-2<<","<<Nz-2<<") = "<<newintgrid.get(Nx-2,Ny-2,Nz-2)<<endl;
		//Find out the indices of each of the positions, assign and print:
		n = 0; vector<Eigen::Vector3i> allindx;
		ofstream ofile; ofile.open(("indices."+to_string(Npos)+".dat"));
		cout<<"Performing INDEXING of 3D double positions:"<<endl;
		ofile<<"# r(x)"<<'\t'<<"r(y)"<<'\t'<<"r(z)"<<'\t'<<"index(i)"<<'\t'<<"index(j)"<<'\t'<<"index(k)"<<endl;
		while(n<Npos){
			Eigen::Vector3i rindx;
			Eigen::Vector3d r; r = allpos[n];
			rindx = position2cell<Nx,Ny,Nz>(r, r_o, d);
			allindx.push_back(rindx);
			ofile<<r(0)<<'\t'<<r(1)<<'\t'<<r(2)<<'\t'<<rindx(0)<<'\t'<<rindx(1)<<'\t'<<rindx(2)<<endl;
			//system ("clear");
			cout<<'\r'<<"Step Number: "<<n<<std::flush;
			//getch();
			n = n + 1;
		}
		ofile.close();
		cout<<endl;
		//move on to the next iteration
		Npos = Npos*2;
		cout<<"|############################################################################|"<<endl;
		cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"<<endl;
	}
	return 0;
}