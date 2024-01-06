#!/usr/bin/python3
from math import sqrt

CRANGE=((30,38),(38,46),(46,54))
NUMFORM=".3f"
\
def magnitude(v): return sqrt(sum([x*x for x in v]))
def parmove(args):
    myf=open(args["infile"],"r")
    XC=float(args["xc"])
    YC=float(args["yc"])
    ZC=float(args["zc"])
    if args["outfile"] is not None: outf=open(args["outfile"],"w")
    else: outf=None
    for l in myf:
        l=l[0:len(l)-1] #Remove "\n" from the end
        if l[0:4]=="ATOM" or l[0:6]=="HETATM":
            foreline=l[0:CRANGE[0][0]]
            xval=float(l[CRANGE[0][0]:CRANGE[0][1]])+XC
            yval=float(l[CRANGE[1][0]:CRANGE[1][1]])+YC
            zval=float(l[CRANGE[2][0]:CRANGE[2][1]])+ZC
            backline=l[CRANGE[2][1]:]
            outline=foreline+format(format(xval,NUMFORM),">8s")+format(format(yval,NUMFORM),">8s")+format(format(zval,NUMFORM),">8s")+backline
            if outf: outf.write(outline+"\n")
            else: print(outline)
        else:
            if outf: outf.write(l+"\n")
            else: print(l)
