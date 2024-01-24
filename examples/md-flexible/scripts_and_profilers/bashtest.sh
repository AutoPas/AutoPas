#!/bin/bash
#SBATCH -J bashtest
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=cm2   #not entirely sure what these 2 lines do yet ^^
#SBATCH --partition=cm2_std
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=200mb
#SBATCH --mail-user=johannes.riemenschneider@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:05:00

module load slurm_setup

echo ./md-flexible --yaml-filename twoMoleculesTrueMultisiteBundled.yaml
./md-flexible --yaml-filename twoMoleculesTrueMultisiteBundled.yaml

