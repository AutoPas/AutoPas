#!/bin/bash
source $(cd "$(dirname "$0")"; pwd -P)/ci_funcs.sh

# load MPI environment
load_MPI

# build
mkdir -p build && cd build

if [[ $? == 0 ]]; then

    #export OMPI_FC=gfortran-8 
    #export OMPI_CC=gcc-8 
    #export OMPI_CXX=g++-8 
    #CC=/usr/lib64/mpi/gcc/openmpi3/bin/mpicc CXX=/usr/lib64/mpi/gcc/openmpi3/bin/mpicxx FC=/usr/lib64/mpi/gcc/openmpi3/bin/mpifort ${CMAKE} .. -DCM_ALL_FORTRAN=ON -DCM_ALL_USE_F08=ON
    ${CMAKE} .. -DCM_ALL_FORTRAN=ON -DCM_ALL_USE_F08=ON
    make VERBOSE=1
fi
