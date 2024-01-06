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
    l=l[0:len(l)-1] #Remove the '\n' in the end
    bpart=str.split(l[0:6],'-')
    p1=str.upper(bpart[0].strip())
    p2=str.upper(bpart[1].strip())
    if (p1 in atoms) and (p2 in atoms):
        kb=float(l[7:13].strip())*1000
        l=float(l[16:23].strip())/10
        print(p1,p2,round(l,5),round(kb,4))
