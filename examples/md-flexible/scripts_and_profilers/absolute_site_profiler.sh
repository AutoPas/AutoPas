#!/bin/bash
#SBATCH -J bashtest
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=200mb
#SBATCH --mail-user=johannes.riemenschneider@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:05:00

#module load slurm_setup

container_prefix="LC" #"VL", "LC"
data_layout="SoA" #"SoA","AoS"
site_counts="1 2 5"
densities="0.5 0.75 1"
functor="BundlingApproach" #"NewFunc","OldFunc", "SingleSiteEmulator", "BundlingApproach"
samples="0 1"

path_to_input_file="./"

for site_count in $site_counts; do
  for density in $densities; do
    for sample in $samples; do
      filename="${container_prefix}${data_layout}Sites${site_count}Density${density}${functor}.yaml"
      echo "Handling ${filename}"
      ./md-flexible --yaml-filename "${filename}"
      echo ""
    done
    echo "-----------------------------------------------"
    echo ""
  done
done