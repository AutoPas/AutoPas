#!/bin/bash

#SBATCH -J Spinodal_Checkpoint
#SBATCH -o ./%x%j.%N.out
#SBATCH -D .
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny

export OMP_NUM_THREADS=28

./build_mp/examples/md-flexible/md-flexible --yaml-filename ./build_mp/examples/md-flexible/SpinodalDecomposition_equilibration.yaml