#ifndef INCLUDED_PHYSICS
#define INCLUDED_PHYSICS 1
namespace physics
{
  void gradientDescentPhysics(const System& sys,Molecule* target,int nsteps=3500,double D=1e-4,double R=2.5e-5,double L=2.25e-5,double L2=5e-5);
  void atomisticGradientDescentPhysics(const System& sys,Molecule* target,int nsteps=10000,double D=3e-8,double DU=1e-6,double DS=1e-7,double R=3e-8,double L=1e-5);
  void optimizeMolecule(Molecule* mol,const ForceField& ff,int nsteps=2500,double D=2.5e-8);
}
#endif
