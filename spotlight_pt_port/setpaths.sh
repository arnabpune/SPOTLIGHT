#!/bin/bash

LOC="$1"
for f in `ls *.cpp`
do
	sed -i '/std::string PTPATH=/c\static std::string PTPATH="'$LOC'";' $f
done
