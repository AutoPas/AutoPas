#!/bin/bash

#SBATCH -J Alpha_Gamma
#SBATCH -o ./%x%j.%N.out
#SBATCH -D .
#SBATCH --clusters=serial
#SBATCH --partition=serial_std

value_array_gamma='0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'
value_array_alpha='0.7 0.75'

export OMP_NUM_THREADS=1

for a in ${value_array_alpha}
do
    (
    for g in ${value_array_gamma}
    do
        (
            export ALPHA=${a}
            export GAMMA=${g}
            strategy="ReinforcementLearning"
            simulation="fallingDrop"
            yaml_file="./build/examples/md-flexible/${simulation}ReinforcementLearning.yaml"
            txt_file="results/alpha_${a}_gamma_${g}.txt"
            ./build/examples/md-flexible/md-flexible "--yaml-filename" ${yaml_file} | tee ${txt_file}
            unset ALPHA
            unset GAMMA
        )
    done
    )
done