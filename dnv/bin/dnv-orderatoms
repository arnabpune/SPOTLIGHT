#!/bin/bash

template=$1
inf=$2

anames=`grep ^ATOM $template|cut -c13-16`
grep "^ *REMARK" $inf
for x in $anames
do
	grep -E "^[HA][ET][TO][AM][T ][M ]...... *$x[ A-Za-z]" $inf|head -1
done
