#/bin/bash

Mols=(32 64 128 256 512 1024 2048 4096 8192)
Reps=(100000 10000 10000 1000 1000 1000 100 10 5)

# iterate over containers
for iCont in {0..1} ;
do
	echo Container $iCont:
	# iterate over molecules with the correct repetition
	for i in {0..8}; do ./md-main $iCont ${Mols[$i]} ${Reps[$i]}; done > output-$iCont.txt
done
