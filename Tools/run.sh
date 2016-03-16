#!/bin/bash

# Check if directory gauss exists, make if not
if [ ! -d "gauss" ]; then
	echo "Directory gauss does not exist. Creating..."
	mkdir gauss
else
	echo "Clearing Gauss..."
	touch gauss/hey
	rm gauss/*
fi

if [ ! -d "IO" ]; then
	echo "Directory IO does not exist. Creating..."
	touch IO/hey
	mkdir IO
else
	echo "Clearing IO..."
	rm IO/*
fi

cp box.txt IO/box.0.in

./MonteCarloOpenMP config.ini 5000000 0 0

cp IO/box.0.out IO/box.1.in

./MonteCarloOpenMP config.ini 0 1000 1


