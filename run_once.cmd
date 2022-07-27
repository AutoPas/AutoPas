#!/bin/bash

#SBATCH -J Alpha_Gamma
#SBATCH -o ./%x%j.%N.out
#SBATCH -D .
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --get-user-env

a='0'
g='0'
export "ALPHA=${a}"
export "GAMMA=${g}"
strategy="ReinforcementLearning"
simulation="fallingDrop"
yaml_file="./build/examples/md-flexible/${simulation}ReinforcementLearning.yaml"
txt_file="alpha_${a}_gamma_${g}.txt"
./build/examples/md-flexible/md-flexible "--yaml-filename" ${yaml_file}
unset "ALPHA"
unset "GAMMA"