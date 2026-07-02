#!/bin/bash
source $(cd "$(dirname "$0")"; pwd -P)/ci_funcs.sh

# load MPI environment
load_MPI

# build
mkdir -p build && cd build

if [[ $? == 0 ]]; then

    #CC=/usr/lib64/mpi/gcc/openmpi3/bin/mpicc CXX=/usr/lib64/mpi/gcc/openmpi3/bin/mpicxx FC=/usr/lib64/mpi/gcc/openmpi3/bin/mpif90 ${CMAKE} ..
    ${CMAKE} .. -DCM_ALL_FORTRAN=ON -DCM_ALL_VORONOI=ON -DCM_ALL_VTK_OUTPUT=ON -DCM_ALL_TESTS=ON -DCM_ALL_AUTO_DOC=ON -DVTK_DIR=/usr/local/lib/cmake/vtk-7.1 
    make VERBOSE=1
    make test
fi
