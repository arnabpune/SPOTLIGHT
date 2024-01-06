#include "./tensor.hpp"

using namespace std;

/*
We will define the a bounding box with an origin and box vectors
*/
pair<Eigen::Vector4d,Eigen::Vector3i> define_bounding_box(vector<Eigen::Vector3d> pos, double rcut){
	//we need to find the minimum and maximum position vectors, and also the d-spacing
	Eigen::Vector3d rmin,rmax;
	double dbox = findlatticespacing(pos, rcut);
	Eigen::Vector3d dr(dbox,dbox,dbox);
	//bearing in mind that the minimum position shifted by the d-vector is the origin
	rmin = (findminpos(pos)) - dr; rmax = (findmaxpos(pos)) + dr;
	cout<<"r_min = "<<rmin.transpose()<<", r_max = "<<rmax.transpose()<<endl;
	Eigen::Vector3d rbox; rbox = rmax-rmin;
	//we then define our box vectors a,b,c to be the rounded off figures for each cell axis
	int a,b,c;
	a = (int)((round(rbox(0)))); b = (int)((round(rbox(1)))); c = (int)((round(rbox(2))));
	//we then create the final datatypes and create the returnable pair and then return.
	Eigen::Vector4d box_o_d(rmin(0),rmin(1),rmin(2),dbox);
	Eigen::Vector3i boxvec(a,b,c);
	pair<Eigen::Vector4d,Eigen::Vector3i> p; p = make_pair(box_o_d,boxvec);
	return p;
}

/*
Assign a position to a cell in a 3D grid:
In the previous function we defined our origin r_o, our box d-spacing d and the box vectors a,b,c.
In this section we are now assuming Nx = a/d, Ny = b/d and Nz = c/d. This allows us to assign a 
unique box to a coordinate in such a PERIODIC box, by using the formula
{i,j,k} = round(r_i-r_o_i)/d,round(r_j-r_o_j)/d,round(r_k-r_o_k)/d  with correction for periodic conditions.
This formula implies that the rounded off distance divided by the d-spacing is equivalent to an index,
or is atleast the nearest index.
The other method would be to search for the nearest point but this is much more complex in timescale
*/

template<int Nx,int Ny,int Nz>
Eigen::Vector3i position2cell(Eigen::Vector3d r, Eigen::Vector3d r_o, double d){
	//first we find the difference
	Eigen::Vector3d rdiff; rdiff = r-r_o;
	//then we use the predefined pbc correction to measure and correct the calculated positions
	Eigen::Vector3i ri; ri = pbc_index_correction<Nx,Ny,Nz>( ((int)(round((rdiff(0))/d))),0,0 );
	Eigen::Vector3i rj; rj = pbc_index_correction<Nx,Ny,Nz>( 0,((int)(round((rdiff(1))/d))),0 );
	Eigen::Vector3i rk; rk = pbc_index_correction<Nx,Ny,Nz>( 0,0,((int)(round((rdiff(2))/d))) );
	//we then extract the values to a meaningful box index and return.
	Eigen::Vector3i boxindx( ri(0), rj(1), rk(2) );
	return boxindx;
}