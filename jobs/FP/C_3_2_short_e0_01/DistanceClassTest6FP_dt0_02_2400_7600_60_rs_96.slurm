#!/bin/bash
#SBATCH -J 6fprs96
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
#SBATCH --time=00:45:00

module load slurm_setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
MYPROG="./md-flexible --yaml-filename DistanceClassTest6FP_dt0_02_2400_7600_60_c_3_2_e0_01_rs_96.yaml"

$MYPROG