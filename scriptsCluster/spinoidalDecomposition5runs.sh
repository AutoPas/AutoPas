#!/bin/bash
#SBATCH -J virialSpin
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ../build
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=112
#SBATCH --export=NONE
#SBATCH --time=08:00:00

module load slurm_setup
module load gcc/15.2.0
module load xerces-c/3.2.1

echo "Running normal..."
    ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/SpinodalDecomposition_equilibration.yaml 

echo "Running normal..."
    ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/SpinodalDecomposition_equilibration.yaml

echo "Running normal..."
    ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/SpinodalDecomposition_equilibration.yaml

echo "Running perf profiling..."
    perf record -g ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/SpinodalDecomposition_equilibration.yaml

echo "Running gprof profiling..."

    gprof ./examples/md-flexible/md-flexible --yaml-filename ./examples/md-flexible/SpinodalDecomposition_equilibration.yaml ./gmon.out > profile_report_gprof.txt
    echo "gprof profiling completed. Output saved to profile_report_gprof.txt."