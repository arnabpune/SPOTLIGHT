#!/bin/bash

x=`ls *_out.pdbqt`
for f in $x
do
	if test -f $f"_mod1.pdbqt"
	then
		echo "Model 1 extracted already. Skipping"
	else
		awk -e '{print; if($0=="ENDMDL") exit}' $f > $f"_mod1.pdbqt"
	fi
	if test -f $f"_red.pdb"
	then
		echo "Already converted. Skipping"
	else
		obabel -h -i pdbqt $f"_mod1.pdbqt" -o pdb -O $f"_red.pdb"
	fi
done
