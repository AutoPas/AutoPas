#!/bin/bash
source $(cd "$(dirname "$0")"; pwd -P)/ci_funcs.sh

# load MPI environment
load_MPI

# build
mkdir -p build && cd build

if [[ $? == 0 ]]; then

    #CC=/usr/lib64/mpi/gcc/openmpi3/bin/mpicc CXX=/usr/lib64/mpi/gcc/openmpi3/bin/mpicxx FC=/usr/lib64/mpi/gcc/openmpi3/bin/mpif90 ${CMAKE} ..
    ${CMAKE} ..
    make VERBOSE=1
fi
