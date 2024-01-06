#!/usr/bin/python3

import sys,os

parms=[]
fn=sys.argv[1]
atf=open(fn,"r")
for l in atf:
    atn=l[0:2].strip()
    wt=l[3:8]
    parms.append((atn,wt))

atf.close()
fn=sys.argv[2]
nbf=open(fn,"r")
for l in nbf:
    l=l[0:len(l)-1] #Remove the '\n' in the end
    atn=l[0:7].strip()
    wt="<weight>"
    for el in parms:
        if el[0]==atn:
            wt=str(el[1])
            break
    sig=round(float(l[12:21].strip())/1.12246,5)
    eps=l[22:].strip()
    print(str.upper(atn),"<valence>",wt,"<charge>",sig,eps,sig/2,"<hybrid>")
