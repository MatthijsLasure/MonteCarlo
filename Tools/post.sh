#!/bin/bash

cat IO/out.0.txt > total.txt

for i in $(seq 1 $1)
do
	sed '1d' "IO/out.$i.txt" >> total.txt
done
