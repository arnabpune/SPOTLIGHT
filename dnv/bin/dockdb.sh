#!/bin/bash

x="$1"
prot=$2
cx=$3
cy=$4
cz=$5
sx=$6
sy=$7
sz=$8
hyd="$9"
if test -f $x".pdb"
then
	echo
else
	wget "https://www.drugbank.ca/structures/small_molecule_drugs/$x.smiles"
	#sed s/"+0$"//g $x".pdb" > temp.pdb
	#mv temp.pdb $x".pdb"
fi
if test -f $x".pdbqt"
then
	echo "pdbqt exists"
else
	cat $x".smiles"|sed s/"@"/""/g|obabel -ismi -o pdb -O $x".pdb" --gen3d
	if test -n "$hyd"
	then
		obabel -ipdb $x".pdb" -opdbqt -O $x".pdbqt" -r -xh
	else
		obabel -ipdb $x".pdb" -opdbqt -O $x".pdbqt" -r
	fi
fi
if test -f $x"_out.pdbqt"
then
	echo "Docking result exists: "$x"_out.pdbqt"
else
	vina --receptor $prot --ligand $x".pdbqt" --center_x $cx --center_y $cy --center_z $cz --size_x $sx --size_y $sy --size_z $sz
	if test $? -ne 0
	then
		echo "Docking failed!! Drug: $x.pdbqt"
		rm $x".pdb" $x".pdbqt"
		read -p "Press enter to continue..." ent
	fi
fi
