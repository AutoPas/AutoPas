#!/bin/bash
#SBATCH -J 1a1
#SBATCH -D ./
#SBATCH -o ./%x.%j.%N.out
#SBATCH -e ./%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --cpus-per-task=112
#SBATCH --time=04:00:00

module load slurm_setup
export OMP_NUM_THREADS=56
export OMP_PLACES="{0:56}"
export OMP_PROC_BIND=close
MYPROG="./md-flexible --yaml-filename LUT_6FPF_dt0_015_60_e1_a1_therm3_5.yaml"

$MYPROG