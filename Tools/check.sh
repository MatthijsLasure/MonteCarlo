#!/bin/bash

while :
do
	line=$(tail -n -1 $1)
	printf "$line \\r"
	sleep 1
done
