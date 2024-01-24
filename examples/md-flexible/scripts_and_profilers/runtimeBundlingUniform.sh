#!/bin/bash
#SBATCH -J runtimeBundling
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --clusters=cm2  #not entirely sure what these 2 lines do yet ^^
#SBATCH --partition=cm2_std
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --mem=200mb
#SBATCH --mail-user=johannes.riemenschneider@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:05:00

#module load slurm_setup

containerIndices=(0 1)
containerTypes=("LinkedCells" "VerletLists" "VerletClusterLists")
traversalOfContainer=("lc_c08" "vl_list_iteration" "vcl_c06")
containerPrefixes=("LC" "VL" "VCL") #"VL", "LC", "VCL"

data_layout="SoA" #"SoA","AoS"
site_counts=(1 2 5)
densities=(0.5 0.75)
num_threads=("4" "8")
functor="BundlingApproach" #"NewFunc","OldFunc", "SingleSiteEmulator", "BundlingApproach"
functorDescriptionInFile="Lennard-Jones Bundling-Approach"
samples="0 1"
newton3States=("disabled" "enabled")

path_to_input_file="./"

createInputFile () {
  #containerIndex data_layout site_count density num_threads functor functorDescriptionInFile

  #build helper strings
  siteTypes="0"
  possibleOtherSitePositions=("filler" "[0.01, 0, 0]" "[-0.01, 0, 0]" "[0, 0.01, 0]" "[0, -0.01, 0]") #filler since index starts at 1
  sitePositions="[0, 0, 0]"

  for ((i=1; i<site_count; i++)); do
      siteTypes+=", 0"
      sitePositions+=", "
      sitePositions+=${possibleOtherSitePositions[i]}
  done
  spacing=$(echo "scale=10; 1/${density}" | bc)
  float_temp_num_particles=$(echo "8*8*8*${density}" | bc)
  temp_number_of_particles=$(printf "%.0f" "$float_temp_num_particles")
  number_of_particles=$(printf "%d" "$temp_number_of_particles")

  echo "# This yaml file is for single-site molecular simulation. Uncomment the Molecules option to run this experiment using
# md-flexible compiled for multi-site molecules.
# container                        :  [LinkedCells, VarVerletListsAsBuild, VerletClusterLists, VerletLists, VerletListsCells, PairwiseVerletLists]
container						 : [${containerTypes[containerIndex]}]
#traversal						 : [${traversalOfContainer[containerIndex]}]
traversal                            :  [ds_sequential, lc_sliced, lc_sliced_balanced, lc_sliced_c02, lc_c01, lc_c01_combined_SoA, lc_c04, lc_c04_HCP, lc_c04_combined_SoA, lc_c08, lc_c18, vcl_cluster_iteration, vcl_c06, vcl_c01_balanced, vcl_sliced, vcl_sliced_balanced, vcl_sliced_c02, vl_list_iteration, vlc_c01, vlc_c18, vlc_sliced, vlc_sliced_balanced, vlc_sliced_c02, vvl_as_built, vlp_c01, vlp_c18, vlp_sliced, vlp_sliced_balanced, vlp_sliced_c02]
data-layout                      :  [${data_layout}] #[AoS, SoA]
# newton3                          :  [disabled, enabled]
newton3                          : [${newton3State}]
verlet-rebuild-frequency         :  10
verlet-skin-radius-per-timestep  :  0.02
verlet-cluster-size              :  4
selector-strategy                :  Fastest-Mean-Value
tuning-strategies                :  [predictive-tuning]
tuning-interval                  :  6000
tuning-samples                   :  10
functor                          :  ${functorDescriptionInFile}
cutoff                           :  2
box-min                          :  [0, 0, 0]
box-max                          :  [9, 9, 9]
cell-size                        :  [1]
deltaT                           :  0.0
iterations                       :  12000
boundary-type                    : [periodic, periodic, periodic]
Sites:
  0:
    epsilon                      :  1.
    sigma                        :  1.
    mass                         :  1.
# Uncomment below to run a multi-site simulation.
Molecules:
  0:
    site-types                   :  [ ${siteTypes} ]
    relative-site-positions      :  [ ${sitePositions} ]
    moment-of-inertia            :  [ 1., 1., 1. ]
Objects:
  CubeUniform:
    0:
      numberOfParticles            :  ${number_of_particles}
      particle-type-id           :  0
      box-length                 :  [8, 8, 8]
      bottomLeftCorner           :  [0.5, 0.5, 0.5]
      velocity                   :  [0, 0, 0]
no-flops                         :  false
no-end-config                    :  true
no-progress-bar                  :  true
vtk-filename                     :  ${filename}
vtk-write-frequency              :  100
" >| "${containerPrefixes[$containerIndex]}${data_layout}Sites${site_count}Density${density}${functor}UniformN3${newton3State}.yaml"
}

for containerIndex in "${containerIndices[@]}"; do
  for site_count in "${site_counts[@]}"; do
    for density in "${densities[@]}"; do
      for threads in "${num_threads[@]}";do
        for newton3State in "${newton3States[@]}";do
          filename="${containerPrefixes[containerIndex]}${data_layout}Sites${site_count}Density${density}${functor}UniformN3${newton3State}.yaml"
          #echo "test" >| notYetThere${containerIndex}_${site_count}_${density}.txt
          createInputFile
          for sample in $samples; do
            echo ""

            echo "Handling ${filename} with ${threads} threads"
            OMP_NUM_THREADS=${threads} ./md-flexible --yaml-filename "${filename}"
            echo ""
          done
          echo "-----------------------------------------------"
          echo ""
        done
      done
    done
  done
done