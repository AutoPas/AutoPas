#!/bin/sh
alpha=1
gamma=1
value_array='0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1'

for a in ${value_array}
do
    (
    for g in ${value_array}
    do
      (
        cd src/autopas/selectors/tuningStrategy || exit
        sed -ie "s/double _alpha = 1/double _alpha = ${a}" ReinforcementLearning.h
        sed -ie "s/double _gamma = 1/double _gamma = ${g}" ReinforcementLearning.h
        cd ../../../../
        cd build-clang
        make autopas
        cd examples/md-flexiable || exit
        mkdir "alpha_${a}_gamma_${g}"
        cd "alpha_${a}_gamma_${g}"


        for simulation in explodingLiquid fallingDrop
        do
          (
          mkdir ${simulation}
          cd ${simulation} || exit
          for strategy in ReinforcementLearning
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
        cd ../../../
        cd src/autopas/selectors/tuningStrategy || exit
        sed -ie "s/double _alpha = ${a}/double _alpha = 1" ReinforcementLearning.h
        sed -ie "s/double _gamma = ${g}/double _gamma = 1" ReinforcementLearning.h
        cd ../../../../
      )
    done
    )
done


