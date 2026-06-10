#!/bin/bash

# $1: Programname
# $2: Outputdirectory (should be relative, for symlnks to be relative)

set -x

PROGRAM=${1##*/}
PROGRAMPATH=${1%/*}

for NP in 4 16 21 64 127 128
do
  #rm -f $2/${PROGRAM}_${NP}.dat
  rm -f $2/${PROGRAM}_${NP}.bin
	mpirun --oversubscribe -np $NP $1 $2 > /dev/null
	#mpirun --oversubscribe -np $NP $1 $2 > $2/${PROGRAM}_${NP}.dat
	#rm -f $2/${PROGRAM}_f_${NP}.dat
	#rm -f $2/${PROGRAM}_f_${NP}.bin
	#ln -s $2/${PROGRAM}_${NP}.dat $2/${PROGRAM}_f_${NP}.dat
	#ln -s $2/${PROGRAM}_${NP}.bin $2/${PROGRAM}_f_${NP}.bin
done
