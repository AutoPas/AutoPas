#!/bin/bash

#SBATCH -J Alpha_Gamma_once
#SBATCH -o ./%x%j.%N.out
#SBATCH -e ./%x%j.%N_error.out
#SBATCH -D .
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --get-user-env
#SBATCH --export=ALPHA,Gamma,OMP_NUM_THREADS

a='0'
g='0'
export OMP_NUM_THREADS=28
export ALPHA=${a}
export GAMMA=${g}
strategy="FullSearch"
simulation="fallingDrop"
yaml_file="./build_mp/examples/md-flexible/${simulation}${strategy}NoProgress.yaml"
txt_file="alpha_${a}_gamma_${g}.txt"
./build_mp/examples/md-flexible/md-flexible "--yaml-filename" ${yaml_file}
unset ALPHA
unset GAMMA