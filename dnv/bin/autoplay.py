import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj

def fileen(fnam):
    f=open(fnam,'r')
    for x in f:
        ens=str.split(str.split(x,' ')[1],'\t')
        en=float(ens[0])
        hen=float(ens[1])
        break
    return en,hen
def filekey(fnam):
    e1,e2=fileen(fnam)
    return e1
#os.chdir("/home/venkata/CADD/betacatenin/results/size19")
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")] #Select all PDB files in this folder
file_names=sorted(file_names,key=filekey,reverse=False)
#rc("open protein.gro")
i=0
rc("2dlabel create energy xpos 0.05 ypos 0.95")
rc("2dlabel create henergy xpos 0.05 ypos 0.9")
ind=0
rc("surf")
while ind<len(file_names):
        fn=file_names[ind]
        print("Loading:",fn)
        en,hen=fileen(fn)
        rc("2dlabel change energy text \"Energy: "+str(en)+"\"")
        rc("2dlabel change henergy text \"Hotspot Energy: "+str(hen)+"\"")
        replyobj.status("Processing " + fn)  #show what file we're working on
        rc("open " + fn)
	#rc("align DNV ~DNV") # put ligand in front of remainder of molecule
        if not i:
            rc("focus: DNV") # center/zoom ligand
            i=1
        rc("surf: DNV") # surface ligand
        #rc("preset apply publication 1") # make everything look nice
        rc("surftransp 55") # make the surface a little bit see-through
        rc("pause")
        rc("close: DNV")
        ind+=1
	# save image to a file that ends in .png rather than .pdb
	#png_name = fn[:-3] + "png"
	#rc("copy file " + png_name + " supersample 3")
	#rc("close all")
