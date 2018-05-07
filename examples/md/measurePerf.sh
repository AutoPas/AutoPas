#!/bin/bash

Mols=(32 64 128 256 512 1024 2048 4096 8192)
Reps=(100000 10000 10000 1000 1000 1000 100 20 20)

# iterate over containers
for iCont in {0..1} ;
do
	echo Container ${iCont}:
	# iterate over molecules with the correct repetition
	echo -e "Number of Molecules\tNumber of Force updates\tElapsed time\tMFUPS\tFLOPs\thit rate\tGFLOP/sec" > output-${iCont}.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]}; done >> output-${iCont}.txt
done
# verlet
iCont=2;
	echo Container ${iCont}:
	# iterate over molecules with the correct repetition
	echo -e "Number of Molecules\tNumber of Force updates\tElapsed time\tMFUPS\tFLOPs\thit rate\tGFLOP/sec" > output-${iCont}-verlet-1-0.0.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 1 0.; done >> output-${iCont}-verlet-1-0.0.txt
    echo -e "Number of Molecules\tNumber of Force updates\tElapsed time\tMFUPS\tFLOPs\thit rate\tGFLOP/sec" > output-${iCont}-verlet-5-0.1.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 5 0.1; done >> output-${iCont}-verlet-5-0.1.txt
    echo -e "Number of Molecules\tNumber of Force updates\tElapsed time\tMFUPS\tFLOPs\thit rate\tGFLOP/sec" > output-${iCont}-verlet-10-0.2.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 10 0.2; done >> output-${iCont}-verlet-10-0.2.txt
    echo -e "Number of Molecules\tNumber of Force updates\tElapsed time\tMFUPS\tFLOPs\thit rate\tGFLOP/sec" > output-${iCont}-verlet-20-0.3.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 20 0.3; done >> output-${iCont}-verlet-20-0.3.txt
