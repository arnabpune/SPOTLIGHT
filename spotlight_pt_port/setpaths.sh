#!/bin/bash

LOC="$1"
for f in `ls *.cpp`
do
	echo "#include <string>" > temp
	echo "std::string PTPATH=\"$LOC\"" >> temp
	sed /"std::string PTPATH="/d $f >> temp
	mv temp $f
done
