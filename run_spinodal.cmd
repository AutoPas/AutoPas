#!/bin/bash

#SBATCH -J Spinodal_Checkpoint
#SBATCH -o ./%x%j.%N.out
#SBATCH -D .
#SBATCH --clusters=serial
#SBATCH --partition=serial_std

./build/examples/md-flexible/md-flexible --yaml-filename SpinodalDecomposition_equilibration.yaml