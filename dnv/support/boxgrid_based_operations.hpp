#include "../core_include.hpp"

using namespace std;

/*
A very important function is to define a periodic coordinate system and to also to define
operations like nearest neighbor searching based on grid points and also to correct a position
based on simple periodic geometery formulae.
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#################Our Coordinate System starts at 0##########################|
|########################And end at (Nmax-1)#################################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/

Eigen::Vector3i pbc_index_correction(int i,int j, int k,int Nx,int Ny,int Nz);
template<int Nx,int Ny,int Nz> inline Eigen::Vector3d pbc_index_correction(int i,int j, int k){return pbc_index_correction(i,j,k,Nx,Ny,Nz);}
Eigen::Vector3i pbc_index_correction(int i,int j, int k,int Nx,int Ny,int Nz){
	int ii = i; int jj = j; int kk = k;
	if(k < 0){
		//k is below 0 so it is negative implying it is k-cells below from Nz = Nz+k
		kk = Nz + k;
		if(j < 0){
			//j is below 0 so it is negative implying it is j-cells below from Ny = Ny+j cell
			jj = Ny + j;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else if(j >= Ny){
			//j is above Ny implying it is at j-Ny cell so that Ny is the 0th cell
			jj = j - Ny;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else{
			//do nothing to j, but i is still up for correction
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
	}
	else if(k >= Nz){
		//k is above Nz so it is positive implying it is in i - Nz cell so that Nz is 0th cell
		kk = k - Nz;
		if(j < 0){
			//j is below 0 so it is negative implying it is j-cells below from Ny = Ny+j cell
			jj = Ny + j;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else if(j >= Ny){
			//j is above Ny implying it is at j-Ny cell so that Ny is the 0th cell
			jj = j - Ny;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else{
			//do nothing to j, but i is still up for correction
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
	}
	else{
		//do nothing to k, but j and i can still be corrected
		if(j < 0){
			//j is below 0 so it is negative implying it is j-cells below from Ny = Ny+j cell
			jj = Ny + j;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else if(j >= Ny){
			//j is above Ny implying it is at j-Ny cell so that Ny is the 0th cell
			jj = j - Ny;
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
		else{
			//do nothing to j, but i is still up for correction
			if(i < 0){
				//i is below 0 so it is negative implying it is i-cells below from Nx = Nx+i cell
				ii = Nx + i;
			}
			else if(i >= Nx){
				//i is above Nx so it is positive implying it is in i-Nx cell so that Nx is 0th cell
				ii = i - Nx;
			}
			else{ii = i;}
		}
	}
	//Now we will contruct the vector3d object
	Eigen::Vector3i v(ii,jj,kk);
	return v;
}
