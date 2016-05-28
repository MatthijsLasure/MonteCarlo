#!/bin/bash

if [ -e $2 ]; then rm $2; fi

sed -e 's/\(.*\)\(!\.*\)/\U\1\E\2/' $1 > temp

while read f;
do
	if echo $f | grep -q \!
	then
		echo "$f" >> $2
	else
		echo "${f^^}" >> $2

	fi

done < temp
