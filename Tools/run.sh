#!/bin/bash
set -e

export OMP_NUM_THREADS=8
export KMP_AFFINITY=scatter

if [ ! -d "IO" ]; then
	echo "Directory IO does not exist. Creating..."
	mkdir IO
fi

if [ ! -d "prep" ]; then
	echo "Directory prep does not exist. Creating..."
	mkdir prep
fi

cp box.txt prep/box.0.in
cp solute.txt prep/solute.0.in
cp solute.txt prep/solute.1.in
cp solute.txt prep/solute.2.in
cp solute.txt prep/solute.3.in
cp solute.txt prep/solute.4.in

# r = 20
./MonteCarloGNU config0.ini 1000000 0 0

./BoxScale 0.9 prep/box.0.out prep/box.1.in

# r = 18
./MonteCarloGNU config0.ini 1000000 0 1

./BoxScale 0.9 prep/box.1.out prep/box.2.in

# r = 16.2
./MonteCarloGNU config0.ini 1000000 0 2

./BoxScale 0.960971428 prep/box.2.out prep/box.3.in

# r = 15.57
./MonteCarloGNU config0.ini 1000000 0 3

# Long Gaussian
cp prep/box.3.out prep/box.4.in
./MonteCarloGNU config0.ini 0 25000 4


# Prep
cp box.txt IO/box.0.out
cp solute.txt IO/solute.0.out

# DOIT
for j in {1..50000};
do
	cp IO/box.$(( j - 1 )).out IO/box.$j.in
	cp IO/solute.$(( j - 1 )).out IO/solute.$j.in

	./MonteCarloGNU config1.ini 0 500 $j
done
