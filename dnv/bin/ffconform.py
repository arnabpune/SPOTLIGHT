#!/usr/bin/python3

#REQUIRES orderatoms.sh
#REQUIRES ffrename.py

import sys
import os

#file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")] #Select all PDB files in this folder
file_names=str.split(sys.argv[1],' ')
if len(sys.argv)>2: ff_file_names=str.split(sys.argv[2],' ')
else: ff_file_names=[]
elen=dict()
K=0
for fl in file_names:
    newfnam=str.split(fl,'.')[0]+"_ffel.pdb"
    if len(ff_file_names)>K: ffnam=ff_file_names[K]
    else: ffnam=str.split(fl,'.')[0]+"_ff.pdb"
    outf=open(newfnam,"w")
    inf=open(fl,"r")
    print(fl,":\t with ff file '"+ffnam+"'")
    K+=1
    for ln in inf:
        ln=ln[0:len(ln)-1] #Remove the '\n' at the end of variable:-ln
        if ln[0:4]=="ATOM" or ln[0:6]=="HETATM":
            elnam=str.strip(ln[len(ln)-4:])
            idv=elen.get(elnam)
            if not idv:
                idv=1
                elen[elnam]=1
            elen[elnam]+=1
            outf.write(ln[0:12]+format(elnam+str(idv),"<4s")+ln[16:]+"\n")
            #outf.write(ln[0:12]+ln[len(ln)-4:]+ln[16:]+"\n")
        else: outf.write(ln+"\n")
    outf.close()
    inf.close()
    os.system("dnv-orderatoms "+ffnam+" "+newfnam+" > temp.pdb")
    os.system("dnv-ffrenumber temp.pdb")
    os.system("mv temp_ffrn.pdb "+newfnam)
    os.system("rm temp.pdb")
    print("Output written to:",newfnam)
