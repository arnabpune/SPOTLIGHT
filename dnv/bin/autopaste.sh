#!/bin/bash

fl=`ls *.pdb`
for f in $fl
do
	y=`head -1 $f`
	e1=`echo $y|cut -d" " -f2`
	e2=`echo $y|cut -d" " -f3`
	if test -f $f"_out.pdbqt"
	then
		ex=`grep "REMARK VINA RESULT" $f"_out.pdbqt"|awk -e '{print $4}'`
		echo $f" "$e1" "$e2" "$ex
	else
		echo "Skipping file: $f. Docking output not found at: $f_out.pdbqt"
	fi
done
