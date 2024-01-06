#!/bin/bash

x=$1 #Folder of old files
y=$2 #Folder of new files

ofiles=`ls $x/*.pdb`
for f in $ofiles
do
	fnam=`basename $f|cut -d"." -f1`
	grep ^CONECT $f >> $y/$fnam"_rewritten.pdb"
	echo Processed $f
done
