#!/bin/bash
#SBATCH -J 6fpf
#SBATCH -D ./
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --time=01:30:00

module load slurm_setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
MYPROG="./md-flexible --yaml-filename DistanceClassTest6FPF_dt0_02_2400_30400_60_c_2_5_2_4_e0_01_perf_ad.yaml"

$MYPROG