#include "draw.h"
#include "../graph/Molecule.hpp"


using draw::Screen;

class MoleculeViewer : public Screen
{
  Molecule* mol=nullptr;
  int BUFFER=5;
  double orX=BUFFER,orY=BUFFER,zoomX=1,zoomY=1,sizeX,sizeY,sizeZ; //Sizes are actual sizes in nanometers
  Eigen::Vector3d cog,cov;
  bool nohyd=true;

public:
  MoleculeViewer() : Screen() {}
  MoleculeViewer(Molecule* m) : MoleculeViewer() {mol=new Molecule(m); centre(mol);}
  MoleculeViewer(int x,int y) : Screen(x,y) {}
  MoleculeViewer(int x,int y,Molecule* m) : Screen(x,y) {mol=new Molecule(m); centre(mol);}

public:
  void centre(Molecule* m)
  {
    cog=m->getCentreOfGeometry();
    for(Atom* a : m->getAtoms()) a->setPosition(a->getPosition()-cog);
    std::vector<double> plane=m->bestPlaneOfFit();
    plane[3]=cog(0)*plane[0]+cog(1)*plane[1]+cog(2)*plane[2];
    Eigen::Vector3d norm(plane[0],plane[1],plane[2]),axis=norm.cross(Eigen::Vector3d(0,0,1)); norm/=norm.norm(); axis/=axis.norm();
    double angle=acos(norm.dot(Eigen::Vector3d(0,0,1)));
    cout << "Angle: "<<angle<<"\n";
    if(abs(angle)>0.05)
    {
      Eigen::AngleAxisd axd(-angle,axis);
      for(Atom* a : m->getAtoms()) a->setPosition(axd*(a->getPosition()));
    }
    std::pair<Eigen::Vector3d,Eigen::Vector3d> bounds=m->getContainer();
    sizeX=(::get<1>(bounds))(0)-(::get<0>(bounds))(0);
    sizeY=(::get<1>(bounds))(1)-(::get<0>(bounds))(1);
    sizeZ=(::get<1>(bounds))(2)-(::get<0>(bounds))(2);
    cout << sizeX<<","<<sizeY<<","<<sizeZ<<"\n";
    zoomX=(columns-2*BUFFER)/sizeX;
    zoomY=(lines-2*BUFFER)/sizeY;
    orX=(::get<0>(bounds))(0);
    orY=(::get<0>(bounds))(1);
    centreScreen();
  }
  void centreScreen()
  {
    double yt=lines/2,xt=columns/2;
    untransform(xt,yt);
    cov=Eigen::Vector3d(xt,yt,0);
    //cog=mol->getCentreOfGeometry();
  }


  void untransform(double& x,double& y)
  {
    x-=BUFFER; y-=BUFFER;
    x/=zoomX; y/=zoomY;
    x+=orX; y+=orY;
  }
  void transform(double& x,double& y)
  {
    x-=orX; y-=orY;
    x*=zoomX; y*=zoomY;
    x+=BUFFER; y+=BUFFER;
  }
  void drawSelf()
  {
    double x1,y1,x2,y2;
    for(Atom* a : mol->getAtoms())
    {
      if(nohyd && a->isHydrogen()) continue;
      std::vector<Atom*> bnds=mol->getBondedAtoms(a);
      x1=a->seek_x(); y1=a->seek_y();
      transform(x1,y1);
      for(Atom* b : bnds)
      {
        if(nohyd && b->isHydrogen()) continue;
        x2=b->seek_x(); y2=b->seek_y();
        transform(x2,y2);
        cout << "("<<x1<<","<<y1<<") -> ("<<x2<<","<<y2<<")\n";
        this->drawLine(x1,y1,x2,y2);
      }
    }
    for(Atom* a : mol->getAtoms())
    {
      if(nohyd && a->isHydrogen()) continue;
      x1=a->seek_x(); y1=a->seek_y();
      transform(x1,y1);
      cout << "("<<x1<<","<<y1<<")\t"<<"'"<<(a->toString()[0])<<"'\n";
      this->plot(x1,y1,a->toString()[0],true);
    }
  }

  //Commands
  void process(std::string fullcommand)
  {
    stringstream ss(fullcommand);
    std::string comm;
    ss >> comm;
    std::string trail=stringfx::drain(ss);
    cout << '\n';
    process(comm,trail);
  }
  void process(std::string command, std::string param)
  {
    stringstream ss(param);
    std::string temp,subcomm;
    double val;
    cout << command << "\t" << param << "\n";
    if(command=="quit" || command=="exit") {exit(0);}
    else if(command=="redraw") {return;}
    else if(command=="zoom")
    {
      ss >> subcomm;
      if(subcomm=="x" || subcomm=="y")
      {
        command+=subcomm;
        param=stringfx::drain(ss);
        process(command,param);
        return;
      }
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) val=2; //Default zoom
      if(subcomm=="in") {if(val<1) val=1/val; zoomX*=val; zoomY*=val;}
      if(subcomm=="out") {if(val<1) val=1/val; zoomX/=val; zoomY/=val;}
    }
    else if(command=="zoomx")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) val=2; //Default zoom
      if(subcomm=="in") {if(val<1) val=1/val; zoomX*=val;}
      if(subcomm=="out") {if(val<1) val=1/val; zoomX/=val;}
    }
    else if(command=="zoomy")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) val=2; //Default zoom
      if(subcomm=="in") {if(val<1) val=1/val; zoomY*=val;}
      if(subcomm=="out") {if(val<1) val=1/val; zoomY/=val;}
    }
    else if(command=="scroll")
    {
      ss >> subcomm;
      if(subcomm=="x" || subcomm=="y")
      {
        command+=subcomm;
        param=stringfx::drain(ss);
        process(command,param);
        return;
      }
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="right";}
      if(subcomm=="right") {orX+=val/zoomX;}
      if(subcomm=="left") {orX-=val/zoomX;}
      if(subcomm=="up") {orY-=val/zoomY;}
      if(subcomm=="down") {orY+=val/zoomY;}
    }
    else if(command=="scrollx")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="left";}
      if(subcomm=="right") {orX+=val/zoomX;}
      if(subcomm=="left") {orX-=val/zoomX;}
    }
    else if(command=="scrolly")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; cout<< "'" << temp << "'\n"; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="up";}
      if(subcomm=="up") {orY-=val/zoomY;}
      if(subcomm=="down") {orY+=val/zoomY;}
    }
    else if(command=="rotatex")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="anticlock";}
      if(subcomm=="clock") {val=-val;}
      if(subcomm=="anticlock") {}
      Eigen::Vector3d axis=Eigen::Vector3d(1,0,0);
      double angle=val*PI/180.0;
      cout << "Angle: "<<angle<<"\n";
      centreScreen();
      if(abs(angle)>0.01)
      {
        Eigen::AngleAxisd axd(angle,axis);
        for(Atom* a : mol->getAtoms()) a->setPosition(axd*(a->getPosition()));
      }
    }
    else if(command=="rotatey")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="anticlock";}
      if(subcomm=="clock") {val=-val;}
      if(subcomm=="anticlock") {}
      Eigen::Vector3d axis=Eigen::Vector3d(0,1,0);
      double angle=val*PI/180.0;
      cout << "Angle: "<<angle<<"\n";
      centreScreen();
      if(abs(angle)>0.01)
      {
        Eigen::AngleAxisd axd(angle,axis);
        for(Atom* a : mol->getAtoms()) a->setPosition(axd*(a->getPosition()));
      }
    }
    else if(command=="rotatez")
    {
      ss >> subcomm;
      if(ss.tellg()==-1) val=0;
      else {ss >> temp; val=std::stod(temp);}
      if(val==0) {val=std::stod(subcomm); subcomm="anticlock";}
      if(subcomm=="clock") {val=-val;}
      if(subcomm=="anticlock") {}
      Eigen::Vector3d axis=Eigen::Vector3d(0,0,1);
      double angle=val*PI/180.0;
      cout << "Angle: "<<angle<<"\n";
      centreScreen();
      if(abs(angle)>0.01)
      {
        Eigen::AngleAxisd axd(angle,axis);
        for(Atom* a : mol->getAtoms()) a->setPosition(axd*(a->getPosition()));
      }
    }
    else if(command=="centre")
    {
      centreScreen();
      double rX=columns/2,rY=lines/2;
      untransform(rX,rY);
      rX-=orX; rY-=orY;
      orX=-rX; orY=-rY;
      //shiftX=-rX; shiftY=cog(1)-rY;
    }
    else if(command=="hidehyd") nohyd=true;
    else if(command=="showhyd") nohyd=false;
  }
};
