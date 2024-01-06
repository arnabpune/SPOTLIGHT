#include "OpenMM.h"
#include <iostream>
#include <fstream>
#include "graph/Molecule.hpp"
using OpenMM::Vec3;

namespace openmm_dnv
{
    const std::string PluginPath("/home/venkata/cpp/openmm/install/lib/plugins/");
    const double EPS_FORCE=1e-3;
    class NoOpenMMContextInitializedException : public std::exception {};

    class OpenMMContextWrapper
    {
        OpenMM::System system;
        OpenMM::Context* cxt=nullptr;
        Molecule* mymol;
        std::vector<OpenMM::Vec3> poses;
    public:
        OpenMMContextWrapper() = delete;
        OpenMMContextWrapper(Molecule* m, const ForceField& ff);

        OpenMM::Context* createContext(const std::string& platform="Reference");
        void mapCoordinatesFromContext();
        void mapCoordinatesToContext();
    };
}

openmm_dnv::OpenMMContextWrapper::OpenMMContextWrapper(Molecule* m, const ForceField& ff)
{
    OpenMM::Platform::loadPluginsFromDirectory(openmm_dnv::PluginPath);
    OpenMM::NonbondedForce* nonbond=new OpenMM::NonbondedForce(); system.addForce(nonbond);
    OpenMM::HarmonicBondForce* bondStretches=new OpenMM::HarmonicBondForce(); system.addForce(bondStretches);
    OpenMM::HarmonicAngleForce* bondBends = new OpenMM::HarmonicAngleForce(); system.addForce(bondBends);
    OpenMM::PeriodicTorsionForce* bondTorsions = new OpenMM::PeriodicTorsionForce(); system.addForce(bondTorsions);

    mymol=m;
    Topology mtop(m,ff);

    //Storing positions
    for(Atom* at : m->getAtoms())
    {
        Eigen::Vector3d pos=at->getPosition();
        poses.push_back(OpenMM::Vec3(pos(0),pos(1),pos(2)));
        system.addParticle(at->seek_mass());
        nonbond->addParticle(at->seek_charge(),at->seek_sigma(),at->seek_epsilon());
        std::cout << at->toString()<<" has mass "<<at->seek_mass()<<" at " << pos.transpose() <<"\n";
    }

    //Adding bonds
    for(int i=0;i<mtop.getNumBonds();i++)
    {
        auto bd=mtop.getBondParametersAt(i);
        bondStretches->addBond(get<0>(get<0>(bd)),get<1>(get<0>(bd)),get<0>(get<1>(bd)),get<1>(get<1>(bd)));
        std::cout << "BOND: "<<get<0>(get<0>(bd))<<"-"<<get<1>(get<0>(bd))<<" with "<<get<0>(get<1>(bd))<<","<<get<1>(get<1>(bd)) <<"\n";
    }

    //Adding angles
    for(int i=0;i<mtop.getNumAngles();i++)
    {
        auto bd=mtop.getAngleParametersAt(i);
        auto atoms=get<0>(bd);
        auto pars=get<1>(bd);
        bondBends->addAngle(atoms[0],atoms[1],atoms[2],pars[0],pars[1]+EPS_FORCE); //Urey bradley omitted
        std::cout << "ANGLE "<<atoms[0]<<"-"<<atoms[1]<<"-"<<atoms[2]<<" with "<<pars[0]<<","<<(pars[1]+EPS_FORCE)<<"\n";
    }

    //Adding dihedrals
    for(int i=0;i<mtop.getNumDihedrals();i++)
    {
        auto bd=mtop.getDihedralParametersAt(i);
        auto atoms=get<0>(bd);
        auto pars=get<1>(bd);
        bondTorsions->addTorsion(atoms[0],atoms[1],atoms[2],atoms[3],pars[2],pars[0],pars[1]+EPS_FORCE);
        std::cout << "TORSION "<<atoms[0]<<"-"<<atoms[1]<<"-"<<atoms[2]<<"-"<<atoms[3]<<" with "<<pars[2]<<","<<pars[0]<<","<<(pars[1]+openmm_dnv::EPS_FORCE)<<"\n";
    }

    nonbond->createExceptionsFromBonds(m->generateBondedPairs(),0.5,0.5);
}

OpenMM::Context* openmm_dnv::OpenMMContextWrapper::createContext(const std::string& plat)
{
    OpenMM::VerletIntegrator integrator(0.001); //Timestep in ps (=4 fs)
    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName(plat);
    cxt = new OpenMM::Context(system,integrator,platform);
    cxt->setPositions(poses);
    return cxt;
}

void openmm_dnv::OpenMMContextWrapper::mapCoordinatesFromContext()
{
    const OpenMM::State& state=cxt->getState(OpenMM::State::Positions);
    std::vector<OpenMM::Vec3> poses=state.getPositions();
}
void openmm_dnv::OpenMMContextWrapper::mapCoordinatesToContext()
{
    std::vector<OpenMM::Vec3> posesx;
    for(Atom* at : mymol->getAtoms())
    {
        Eigen::Vector3d pos=at->getPosition();
        posesx.push_back(OpenMM::Vec3(pos(0),pos(1),pos(2)));
    }
    cxt->setPositions(posesx);
}
