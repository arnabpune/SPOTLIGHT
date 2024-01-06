#include <cmath>
#include <iostream>
#include <vector>
using namespace std;
#pragma once

/*const double PI=3.141592653589,beta = 4.17000;//minus_one_by_beta=-1.0/beta;
double throwarandompoint(double min=0, double max=1){
	random_device rd; mt19937 eng(rd());
	uniform_real_distribution<> distr(min,max);
	return distr(eng);
}
double throwarandompoint_normal(double mean=0,double stdev=1.0) //PCMODIFIED: ADDITION (added mean and stdev parameters)
{
	random_device rd; mt19937 eng(rd());
	normal_distribution<double> distr(mean,stdev);
  return distr(eng);
}
inline double boltzmannFactor(const double& energy,const double& temp) {return ::pow(2.71828,(-beta*energy)/100.00);}
template<class T> const T& randomSelect(const std::vector<T>& v) {return v[(int)(throwarandompoint(0,v.size()))];}*/
namespace maths
{
  template<class T> T min(const std::vector<T>& v)
  {
    if(v.size()<=0)
      return T();
    T m=v[0];
    for(int i=1;i<v.size();i++) {if(m>v[i]) m=v[i];}
    return m;
  }
  template<class T> T max(const std::vector<T>& v)
  {
    if(v.size()<=0) return T(0);
    T m=v[0];
    for(int i=1;i<v.size();i++) {if(m<v[i]) m=v[i];}
    return m;
  }
}
