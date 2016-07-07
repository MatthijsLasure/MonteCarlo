#!/bin/bash

while :
do
	file=$( grep '^OUT' $1 | tail -1 | tr -s ' ' | cut -d ' ' -f2)
	line=$(tail -n -1 $file)
	printf "$file $line \\r"
	sleep 1
done
