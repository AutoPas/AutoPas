#!/bin/bash

for i in {0.25,0.50,0.75,1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.00,3.25,3.50}
do
    sed -i "13s/.*/      particle-spacing           :  ${i}/" ./Cubes_hwy.yaml
    echo "Running Cube experiment with ${i} spacing"
    OMP_NUM_THREADS=1 ./md-flexible --yaml-filename Cubes_hwy.yaml --no-progress-bar | tee ./output/benchmark_Cube_${i}.txt
done