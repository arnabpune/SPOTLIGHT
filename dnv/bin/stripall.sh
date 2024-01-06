#!/bin/bash

x=`ls result_prot*.pdb`
for f in $x
do
	fn=`echo $f|cut -d"." -f1`
	if test -f "edit_"$fn"_noH.pdb"
	then
		echo "Converted molecule exists. Skipping: "$fn
	else
		obabel -d -ipdb $f -opdb -O "edit_"$fn"_noH.pdb"
	fi
done
