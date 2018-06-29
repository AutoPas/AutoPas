#!/bin/bash

Mols=(32 64 128 256 512 1024 2048 4096 8192)
Reps=(100000 10000 10000 1000 1000 1000 100 10 5)
cont_to_name=("lc" "direct" "verlet")
func_to_name=("density" "hydro")


for func in {0..1};
do
    echo ----- functor $func
    for iCont in {0..2} ;
    do
        echo Container ${iCont}:
        filename="output-${cont_to_name[$iCont]}-${func_to_name[$func]}.txt"

        echo -e "numparticles\tnumIterations\tMFUPS_AOS\tMFUPS_SOA" > ${filename}
        # iterate over molecules with the correct repetition
        for i in {0..8}; do ./sph-diagram-generation ${Mols[$i]} ${Reps[$i]} ${iCont} ${func} ; done >> ${filename}

    done
done