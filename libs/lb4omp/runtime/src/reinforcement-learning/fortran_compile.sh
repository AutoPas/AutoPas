#!/bin/bash

ml intel/2022a
export LD_LIBRARY_PATH=$(pwd):$LD_LIBRARY_PATH

#RL_files=$(find . -type f -name "*.cpp")
RL_files=$(find . -type f -name "*.cpp" | grep -v "./c_connector.cpp")

echo $RL_files

mpiicpc -lstdc++ -fpic -shared $RL_files -o lib_fortran_rl.so
mpiicpc -lstdc++ -fpic -shared fortran_connector.cpp -L. lib_fortran_rl.so -o lib_fortran_connector.so

