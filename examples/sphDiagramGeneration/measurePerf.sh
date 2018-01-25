#/bin/bash

Mols=(32 64 128 256 512 1024 2048 4096 8192)
Reps=(100000 10000 10000 1000 1000 1000 100 10 5)


for iCont in {0..1} ;
do
	echo Container $iCont:
    echo -e "numparticles\tnumIterations\telapsedTime\tMFUPS\tFLOPs\thit-rate\tGFLOP/s" > output-$iCont.txt
    # iterate over molecules with the correct repetition
    for i in {0..8}; do ./sph-diagram-generation ${Mols[$i]} ${Reps[$i]} $iCont 0 ; done >> output-$iCont.txt

done