#!/bin/bash
#SBATCH -J vtuneBashtest
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
module load intel-oneapi-vtune
export MPS_STAT_LEVEL=4


# Collection
#aps --result-dir=. ./md-flexible --yaml-filename twoMoleculesTrueMultisiteBundled.yaml
#/dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/x86_64/intel-oneapi-vtune/2021.7.1-gcc-ihh6yuf/vtune/2021.7.1/bin64/aps --result-dir=../../../vtune_results ./md-flexible --yaml-filename twoMoleculesTrueM>
#vtune -collect memory-access -result-dir ../../../vtune_results -- ./md-flexible --yaml-filename twoMoleculesTrueMultisiteBundled.yaml
vtune -collect memory-access -data-limit=0 -r ../../../vtune_results -- ./md-flexible --yaml-filename twoMoleculesTrueMultisiteBundled.yaml

# Report
vtune -report summary -r ../../../vtune_results
#aps-report .  # Creates a useful .html report that can be viewed with any browser
#aps-report -a .  # Prints all available stats to stdout
#/dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/x86_64/intel-oneapi-vtune/2021.7.1-gcc-ihh6yuf/vtune/2021.7.1/bin64/aps-report ../../../vtune_results
#/dss/dsshome1/lrz/sys/spack/release/22.2.1/opt/x86_64/intel-oneapi-vtune/2021.7.1-gcc-ihh6yuf/vtune/2021.7.1/bin64/aps-report -a ../../../vtune_results
