#!/bin/bash
#SBATCH -J autopas_highway
#SBATCH -o ./output/%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=end
#SBATCH --mail-user=luis.gall@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:02:00
module load slurm_setup
module load gcc
module load cmake
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./build/examples/md-flexible/md-flexible --yaml-filename ./build/examples/md-flexible/Cube_20_hwy.yaml