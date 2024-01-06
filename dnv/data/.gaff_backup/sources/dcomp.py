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
    bpart=str.split(l[0:12],'-')
    p1=str.upper(bpart[0].strip())
    p2=str.upper(bpart[1].strip())
    p3=str.upper(bpart[2].strip())
    p4=str.upper(bpart[3].strip())
    if (p1=="X" or p2=="X" or p3=="X" or p4=="X") or ((p1 in atoms) and (p2 in atoms) and (p3 in atoms) and (p4 in atoms)):
        m=int(l[13:16].strip())
        ka=float(l[18:25].strip())
        a=float(l[30:39].strip())
        print(p1,p2,p3,p4,round(a,5),round(ka,4),m)
