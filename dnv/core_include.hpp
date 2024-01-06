#ifndef CORE_INCLUDE_HPP
#define CORE_INCLUDE_HPP
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <random>
#include "./Eigen337/Eigen/Core"
#include "./Eigen337/Eigen/Dense"
#include "./Eigen337/Eigen/Geometry"


using namespace std;

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|####################Define a random number generator########################|
|###########################for INT and DOUBLES##############################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/
int genrandomint(int min,int max){
	random_device device;
    mt19937 generator(device());
    uniform_int_distribution<int> distribution(min,max);
	return distribution(generator);
}

double genrandomdouble(double min,double max){
	random_device device;
    mt19937 generator(device());
    uniform_real_distribution<double> distribution(min,max);
	return distribution(generator);
}

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|####################Define a random Position generator######################|
|###########################for INT and DOUBLES##############################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/
Eigen::Vector3i genradnomvec3i(int xm,int xM, int ym, int yM, int zm, int zM){
	Eigen::Vector3i r((genrandomint(xm,xM)), (genrandomint(zm,zM)), (genrandomint(zm,zM)));
	return r;
}

Eigen::Vector3d genrandomvec3d(double xm, double xM, double ym, double yM, double zm, double zM){
	Eigen::Vector3d r((genrandomdouble(xm,xM)), (genrandomdouble(zm,zM)), (genrandomdouble(zm,zM)));
	return r;
}

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#######################Define a min/max finder for a########################|
|#######################vector of 3D-vector positions########################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/

Eigen::Vector3d findminpos(vector<Eigen::Vector3d> poses){
	//consider the absolute zero as the reference point
	Eigen::Vector3d minpos;
	minpos = poses[0];
	int c = 0;//c signifies the component index, 0 = x,1 = y,2 = z
	while(c<3){
		vector<Eigen::Vector3d>::iterator it; it = poses.begin();
		while(it != poses.end()){
			Eigen::Vector3d rc; rc = *it;
			if(minpos(c) > rc(c)){
				minpos(c) = rc(c);//update the cth component
				++it;
			}
			else{ ++it; }
		}
		c = c + 1;
	}
	//return what we got
	return minpos;
}

Eigen::Vector3d findmaxpos(vector<Eigen::Vector3d> poses){
	//consider the absolute zero as the reference point
	Eigen::Vector3d maxpos;
	maxpos = poses[0];
	int c = 0;//c signifies the component index, 0 = x,1 = y,2 = z
	while(c<3){
		vector<Eigen::Vector3d>::iterator it; it = poses.begin();
		while(it != poses.end()){
			Eigen::Vector3d rc; rc = *it;
			if(maxpos(c) < rc(c)){
				maxpos(c) = rc(c);//update the cth component
				++it;
			}
			else{ ++it; }
		}
		c = c + 1;
	}
	//return what we got
	return maxpos;
}

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|############################################################################|
|#######################Define a Monte Carlo algorithm#######################|
|####################to find an average pairwise distance####################|
|############################################################################|
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
*/
/*
We can add a switchingfuncion later
*/
double findlatticespacing(vector<Eigen::Vector3d> pos, double rcut){
	//we will define a random number of trials
	int Ntrials = genrandomint(1,pos.size());
	int n = 0; double sum = 0.0; double nsum = 0.0;
	while(n < Ntrials){
		int indx = genrandomint(0,((pos.size())-1));
		//find it's nearest neighbors
		Eigen::Vector3d r_i; r_i = pos[indx];
		int j = 0;
		while(j < pos.size()){
			Eigen::Vector3d r_j; r_j = pos[j];
			if(j != indx){
				Eigen::Vector3d distvec; distvec = r_i - r_j;
				double d = distvec.norm();
				if(islessequal(d,rcut)!=false){
					//once found sum up the distance of the nearest neighbors
					sum = sum + d; nsum = nsum + 1.0;
					j = j + 1;
				}
				else{j = j + 1;}
			}
			else{j = j + 1;}
		}
		n = n + 1;
	}
	//next divide the distance sum by the number of total trials.
	cout<<"Lattice Spacing Algorithm: SUM = "<<sum<<", NSUM = "<<nsum<<endl;
	double ret = sum/nsum;
	return ret;
}
#endif
