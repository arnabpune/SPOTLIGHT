#ifndef INCLUDED_DRAW
#define INCLUDED_DRAW 1
//#include "commons/commons.h"
#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
namespace draw
{
  inline static char superpose(char c1,char c2) {return c2;}
  class Screen
  {
  protected:
    int lines,columns,buffer=0,gscale=1;
    char FILL='#',EMPTY=' ';
    std::vector<char*> screen;
  public:
    Screen() : Screen(50,50) {}
    Screen(int w,int h)
    {
      columns=w; lines=h;
      screen=std::vector<char*>();
      for(int i=0;i<h;i++) screen.push_back(new char[w]);
      for(int i=0;i<h;i++) {for(int j=0;j<w;j++) screen[i][j]=EMPTY;}
    }
    ~Screen() {for(int i=0;i<lines;i++) delete[] screen[i];}

    inline int getWidth() const {return columns;}
    inline int getColumns() const {return getWidth();}
    inline int getLines() const {return lines;}
    inline int getRows() const {return getLines();}
    inline int getBuffer() const {return buffer;}
    inline void setBuffer(double b) {buffer=b;}

    void draw(std::ostream& os=std::cout)
    {
      predraw(os);
      for(int i=lines-1;i>=0;i--) { os << "\n"; for(int j=0;j<columns;j++) os << screen[i][j];}
      postdraw(os);
      for(int i=0;i<buffer;i++) os << "\n";
    }

  private:
    void set(int x,int y,bool v,bool leak=false)
    {
      if(leak)
      {
        if(x<0 || x>=getRows()) return;
        if(y<0 || y>=getColumns()) return;
      }
      screen[x][y]=(v)?FILL:EMPTY;
    }
    void set(int x,int y,char v,bool leak=false)
    {
      if(leak)
      {
        if(x<0 || x>=getRows()) return;
        if(y<0 || y>=getColumns()) return;
      }
      screen[x][y]=v;
    }
    char get(int x,int y,bool leak=false) const
    {
      if(leak)
      {
        if(x<0 || x>=getRows()) return EMPTY;
        if(y<0 || y>=getColumns()) return EMPTY;
      }
      return screen[x][y];
    }

  public:
    void plot(int y,int x,bool v,bool leak=false) //Variables are badly named. Give (x,y) as input
    {
      if(leak)
      {
        if(x<0 || x>=getRows()) return;
        if(y<0 || y>=getColumns()) return;
      }
      screen[x][y]=(v)?FILL:EMPTY;
    }
    void plot(int y,int x,char v,bool leak=false) //Variables are badly named. Give (x,y) as input
    {
      if(leak)
      {
        if(x<0 || x>=getRows()) return;
        if(y<0 || y>=getColumns()) return;
      }
      screen[x][y]=v;
    }
    char at(int y,int x) const {return screen[x][y];} //Variables are badly named. You will get (x,y) as input

    virtual void predraw(std::ostream& os) {}
    virtual void postdraw(std::ostream& os) {}

    //Drawing methods
    void drawLine(int y1,int x1,int y2,int x2,char chare=' ',bool trace=false) //Variables are badly names (Give input as x1,y1,x2,y2)
    {
      double slope=((double)(x2-x1))/((y2-y1)),cutoff=((double)(lines))/((columns));
      std::cout << "Cutoff: "<<cutoff<<"\tSlope: "<<slope<<"\n";
      char fillc=chare;
      if(chare==EMPTY)
      {
        if(abs(slope)<cutoff) fillc='=';
        else if(abs(slope)>1/cutoff) fillc='|';
        else
        {
          if(slope<1) fillc='/';
          else fillc='\\';
        }
      }
      //slope*=gscale;
      set(x1,y1,superpose(get(x1,y1,true),fillc),true);
      int yp;
      int start=1,end=x2-x1;
      for(int i=start;i<end;i++)
      {
        yp=(int)(i/slope)+y1;
        set(x1+i,yp,superpose(get(x1+i,yp,true),fillc),true);
      }
      start=1; end=y2-y1;
      for(int i=start;i<end;i++)
      {
        yp=(int)(i*slope) +x1;
        set(yp,y1+i,superpose(get(yp,y1+i,true),fillc),true);
      }
    }

    inline void clear() {for(int i=0;i<lines;i++) for(int j=0;j<columns;j++) screen[i][j]=EMPTY;}
  };
}
#endif
