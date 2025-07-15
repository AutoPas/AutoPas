#!/bin/bash
#SBATCH -J 6fprs1
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
#SBATCH --time=00:05:00

module load slurm_setup
export OMP_NUM_THREADS=56
export OMP_PLACES="{0:56}"
export OMP_PROC_BIND=close
MYPROG="./md-flexible --yaml-filename DistanceClassTest6FP_dt0_015_3600_11400_60_c_2_5_2_4_e1_therm3_5_ibi_perf_rs_1.yaml"

$MYPROG