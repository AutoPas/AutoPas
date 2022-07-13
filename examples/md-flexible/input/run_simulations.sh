#!/bin/sh
mkdir "$(date +"%Y-%m-%d_%H-%M")"
cd "$(date +"%Y-%m-%d_%H-%M")"


for simulation in explodingLiquid fallingDrop fallingDrop2 SpinodalDecomposition
do
  (
  mkdir ${simulation}
  cd ${simulation} || exit
  for strategy in FullSearch Predictive ReinforcementLearning
  do
    (
    mkdir ${strategy}
    cd ${strategy} || exit
    for i in 0 1 2
    do
      (
      mkdir "run"${i}
      cd "run"${i} || exit
      yaml_file="../../../../${simulation}${strategy}NoProgress.yaml"
      txt_file="../${strategy}1-${i}.txt"
      ./../../../../md-flexible "--yaml-filename" ${yaml_file} | tee ${txt_file}
      cd ..
      )
    done
    cd ..
    )
  done
  cd ..
  )
done