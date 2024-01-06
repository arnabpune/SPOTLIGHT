#!/bin/bash

x=$1
cap=0
terminal=0
branch=0
hydrogens="HH `grep ^ATOM $x|grep "H *$"|awk -e '{print $2}'` HH"
echo $hydrogens
for l in `grep ^CONECT $x|sed s/"  *"/"_"/g`
do
	an=`echo $l|cut -d"_" -f2`
	echo $hydrogens|grep " $an " 1> /dev/null 2> /dev/null
	if test $? -eq 0
	then
		echo Skipped $an
		continue
	fi
	n=`echo $l|grep -o "_"|wc -l`
	m=$n
	if test $n -eq 2
	then
		cap=`expr $cap + 1`
	else
		i=`expr $n + 1`
		while test $i -gt 2
		do
			an=`echo $l|cut -d"_" -f$i`
			echo $hydrogens|grep " $an " 1> /dev/null 2> /dev/null
			if test $? -eq 0
			then
				n=`expr $n - 1`
			fi
			i=`expr $i - 1`
		done
		if test $n -ge 3
		then
			branch=`expr $branch + 1`
		else
			terminal=`expr $terminal + 1`
		fi
	fi
	echo Processed $an from \"$l\" with n going from n=$m to n=$n 
done
echo
echo "------------------------"
echo
echo "Cap: $cap"
echo "Terminal: $terminal"
echo "Branch (Includes rings): $branch"
