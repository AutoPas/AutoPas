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

containerIndices=(1 0)
containerTypes=("LinkedCells" "VerletLists" "VerletClusterLists")
traversalOfContainer=("lc_c08" "vl_list_iteration" "vcl_c06")
containerPrefixes=("LC" "VL" "VCL") #"VL", "LC", "VCL"
newton3OfContainer=("enabled" "disabled" "disabled")
newton3OfContainerBool=(1 0 1)

data_layout="SoA" #"SoA","AoS"
#site_counts=(1 2 5)
#densities=(0.5 0.75)
num_threads=("4")
functor="OldFunc" #"NewFunc","OldFunc", "SingleSiteEmulator", "BundlingApproach"
functorDescriptionInFile="Lennard-Jones" # Bundling-Approach"
samples="0 1"
#newton3States=("disabled" "enabled")

declare -A densityToMolCount
densityToMolCount["0.1"]=935
densityToMolCount["0.2"]=1715
densityToMolCount["0.3"]=2520
densityToMolCount["0.4"]=3560
densityToMolCount["0.5"]=4311
densityToMolCount["0.6"]=5301
densityToMolCount["0.7"]=5920
densityToMolCount["0.8"]=6970
densityToMolCount["0.9"]=7700
densityToMolCount["1.0"]=8316

path_to_input_file="./"

createInputFile () {
  #containerIndex data_layout site_count density num_threads functor functorDescriptionInFile
  local containerIndexLocal=("${!1}")
  local dataLayoutLocal=("${!2}")
  local siteCountLocal=("${!3}")
  local densityLocal=("${!4}")
  local functorLocal=("${!5}")
  local filenameLocal=("${!6}")
  local filenameWithoutYamlLocal=("${!7}")

  #build helper strings
  siteTypes="0"
  possibleOtherSitePositions=("filler" "[0.01, 0, 0]" "[-0.01, 0, 0]" "[0, 0.01, 0]" "[0, -0.01, 0]" "[0, 0, 0.01]" "[0, 0, -0.01]"
                                    "[0.01, 0.01, 0]" "[0.01, -0.01, 0]" "[-0.01, 0.01, 0]" "[-0.01, -0.01, 0]") #filler since index starts at 1
  sitePositions="[0, 0, 0]"

  for ((i=1; i<site_count; i++)); do
      siteTypes+=", 0"
      sitePositions+=", "
      sitePositions+=${possibleOtherSitePositions[i]}
  done
  #spacing=$(echo "scale=10; 1/${density}" | bc)
  spacing=$(echo "scale=10; (e (l(2)/6))/ (e (l(${density})/3))" | bc -l) #this term is computing 2^(1/6)*density^(-1/3). the exponent (-1/3) was missing in the previous version, that was a mistake

  echo "# This yaml file is for single-site molecular simulation. Uncomment the Molecules option to run this experiment using
# md-flexible compiled for multi-site molecules.
# container                        :  [LinkedCells, VarVerletListsAsBuild, VerletClusterLists, VerletLists, VerletListsCells, PairwiseVerletLists]
container						 : [${containerTypes[containerIndexLocal]}]
#traversal                            :  [ds_sequential, lc_sliced, lc_sliced_balanced, lc_sliced_c02, lc_c01, lc_c01_combined_SoA, lc_c04, lc_c04_HCP, lc_c04_combined_SoA, lc_c08, lc_c18, vcl_cluster_iteration, vcl_c06, vcl_c01_balanced, vcl_sliced, vcl_sliced_balanced, vcl_sliced_c02, vl_list_iteration, vlc_c01, vlc_c18, vlc_sliced, vlc_sliced_balanced, vlc_sliced_c02, vvl_as_built, vlp_c01, vlp_c18, vlp_sliced, vlp_sliced_balanced, vlp_sliced_c02]
traversal						 : [${traversalOfContainer[containerIndexLocal]}]
data-layout                      :  [${dataLayoutLocal}] #[AoS, SoA]
# newton3                          :  [disabled, enabled]
newton3                          : [${newton3OfContainer[containerIndexLocal]}]
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
box-max                          :  [21, 21, 21]
cell-size                        :  [2]
deltaT                           :  0.0
iterations                       :  100
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
  CubeClosestPacked:
    0:
      particle-type-id           :  0
      box-length                 :  [20, 20, 20]
      bottomLeftCorner           :  [0.5, 0.5, 0.5]
      particle-spacing           :  ${spacing}
      velocity                   :  [0, 0, 0]
no-flops                         :  false
no-end-config                    :  true
no-progress-bar                  :  true
vtk-filename                     :  ${filenameWithoutYamlLocal}
vtk-write-frequency              :  100
" >| ${filenameLocal}
}

takeMeasurementsAccordingToParameters(){
  local site_counts_local=("${!1}")
  local densities_local=("${!2}")

  for site_count in "${site_counts_local[@]}"; do
    for density in "${densities_local[@]}"; do
      for threads in "${num_threads[@]}";do
        #for newton3State in "${newton3States[@]}";do
        filename="${containerPrefixes[$containerIndex]}${data_layout}Sites${site_count}TrueDensity${density}${functor}Newton3${newton3OfContainer[$containerIndex]}.yaml"
        filenameWithoutYaml="${containerPrefixes[$containerIndex]}${data_layout}Sites${site_count}TrueDensity${density}${functor}Newton3${newton3OfContainer[$containerIndex]}"
        #echo "test" >| notYetThere${containerIndex}_${site_count}_${density}.txt
        createInputFile containerIndex data_layout site_count density functor filename filenameWithoutYaml
        for sample in $samples; do
          echo ""

          echo "Handling ${filename} with ${threads} threads"
          OMP_NUM_THREADS=${threads} ./md-flexible --yaml-filename "${filename}"
          echo "${site_count},${densityToMolCount[$density]},${threads},${newton3OfContainerBool[containerIndex]}"
          echo ""
        done
        echo "-----------------------------------------------"
        echo ""
        #done
      done
    done
  done
}

echo ""
echo "------------------------------------------------------------"
echo "             Verlet Lists iterating over densities"
echo "------------------------------------------------------------"
echo ""
site_counts=(5)
densities=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
containerIndex=1  #VerletList
takeMeasurementsAccordingToParameters site_counts[@] densities[@]
echo ""
echo "------------------------------------------------------------"
echo "             Verlet Lists iterating over site counts"
echo "------------------------------------------------------------"
echo ""
site_counts=(1 2 3 4 5 6 7 8)
densities=(0.8)
containerIndex=1  #VerletList
takeMeasurementsAccordingToParameters site_counts[@] densities[@]
echo ""
echo "------------------------------------------------------------"
echo "             Linked Cells iterating over densities"
echo "------------------------------------------------------------"
echo ""
site_counts=(5)
densities=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
containerIndex=0  #Linked Cells
takeMeasurementsAccordingToParameters site_counts[@] densities[@]
echo ""
echo "------------------------------------------------------------"
echo "            Linked Cells iterating over site counts"
echo "------------------------------------------------------------"
echo ""
site_counts=(1 2 3 4 5 6 7 8)
densities=(0.8)
containerIndex=0  #Linked Cells
takeMeasurementsAccordingToParameters site_counts[@] densities[@]