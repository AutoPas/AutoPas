#!/bin/bash

ml intel/2022a

RL_files=$(find . -type f -name "*.cpp" | grep -v "./fortran_connector.cpp")

mpiicpc -lstdc++ -fpic -shared $RL_files -o lib_rl.so
mpiicpc -lstdc++ -fpic -shared c_connector.cpp -L. lib_rl.so -o lib_connector.so

export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH