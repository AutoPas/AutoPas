#!/bin/sh
#PBS -l nodes=8

HARMONY_HOME=/usr/local

#PATH=/usr/local/bin:$PATH
#export PATH

HARMONY_S_HOST=brood00
export HARMONY_S_HOST

#HARMONY_S_PORT=1979
#export PATH HARMONY_S_PORT

cd $HARMONY_HOME/example/code_generation
mpirun -np 8 ./gemm /tmp/codegen-gemm TARGET_URL=ssh://brood00//hivehomes/rchen/scratch/code SERVER_URL=ssh://mashie.cs.umd.edu//tmp/codegen-tmp
