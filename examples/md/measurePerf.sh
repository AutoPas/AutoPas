#!/bin/bash

Mols=(32 64 128 256 512 1024 2048 4096 8192)
Reps=(100000 10000 10000 1000 1000 1000 100 20 20)

# iterate over containers
for iCont in {0..1} ;
do
	echo Container ${iCont}:
	# iterate over molecules with the correct repetition
	echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]}; done | tee output-${iCont}.txt
	echo Container ${iCont} soa:

	echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}.txt
    for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} soa; done | tee output-${iCont}-soa.txt

done


# verlet
iCont=2;
	echo "Container ${iCont}:"
	# iterate over molecules with the correct repetition
	echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-1-0.0.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 1 0.; done | tee output-${iCont}-verlet-1-0.0.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-5-0.1.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 5 0.1; done | tee output-${iCont}-verlet-5-0.1.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-10-0.2.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 10 0.2; done | tee output-${iCont}-verlet-10-0.2.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-20-0.3.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 20 0.3; done | tee output-${iCont}-verlet-20-0.3.txt

	echo "Container ${iCont} soa:"
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-1-0.0-soa.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 1 0. soa; done | tee output-${iCont}-verlet-1-0.0-soa.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-5-0.1-soa.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 5 0.1 soa; done | tee output-${iCont}-verlet-5-0.1-soa.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-10-0.2-soa.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 10 0.2 soa; done | tee output-${iCont}-verlet-10-0.2-soa.txt
    echo -e "\"Number of Molecules\"\t\"Number of Force updates\"\t\"Elapsed time\"\t\"MFUPS\"\t\"FLOPs\"\t\"hit rate\"\t\"GFLOP/sec\"" | tee output-${iCont}-verlet-20-0.3-soa.txt
	for i in {0..8}; do ./md-main ${iCont} ${Mols[$i]} ${Reps[$i]} 20 0.3 soa; done | tee output-${iCont}-verlet-20-0.3-soa.txt
