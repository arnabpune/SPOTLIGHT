#!/usr/bin/python3

import sys,os

atoms=[]
fn=sys.argv[1]
atf=open(fn,"r")
for l in atf:
    atn=l[0:2].strip()
    atoms.append(atn)
atf.close()

fn=sys.argv[2]
nbf=open(fn,"r")
for l in nbf:
    #l=l[0:len(l)-1] #Remove the '\n' in the end
    bpart=str.split(l[0:9],'-')
    p1=str.upper(bpart[0].strip())
    p2=str.upper(bpart[1].strip())
    p3=str.upper(bpart[2].strip())
    if (p1 in atoms) and (p2 in atoms) and (p3 in atoms):
        ka=float(l[11:17].strip())
        a=float(l[22:29].strip())
        print(p1,p2,p3,round(a,5),round(ka,4),"0.0","0.0")
