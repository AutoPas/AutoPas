#!/bin/bash

for i in {1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9}
do
    sed -i "13s/.*/      particle-spacing           :  ${i}/" ./Cubes_hwy_A.yaml
    echo "Running Cube experiment with ${i} spacing"
    OMP_NUM_THREADS=1 ./md-flexible --yaml-filename Cubes_hwy_A.yaml --vectorizationPattern 1xVectorLength --no-progress-bar | tee ./output/benchmark_Cube_${i}_1xVec.txt
    OMP_NUM_THREADS=1 ./md-flexible --yaml-filename Cubes_hwy_A.yaml --vectorizationPattern 2xVectorLengthDiv2 --no-progress-bar | tee ./output/benchmark_Cube_${i}_2xVecDiv2.txt
    OMP_NUM_THREADS=1 ./md-flexible --yaml-filename Cubes_hwy_A.yaml --vectorizationPattern VectorLengthDiv2x2 --no-progress-bar | tee ./output/benchmark_Cube_${i}_VecDiv2x2.txt
    OMP_NUM_THREADS=1 ./md-flexible --yaml-filename Cubes_hwy_A.yaml --vectorizationPattern VectorLengthx1 --no-progress-bar | tee ./output/benchmark_Cube_${i}_Vecx1.txt
done