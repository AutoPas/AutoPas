#!/bin/bash
if [ $# -eq 0 ]; then
    moment="$(date '+%y-%m-%d-%k-%M-%S')"
    git clone https://github.com/AutoPas/AutoPas.git -b noah-zonal-triwise "ZonalMD-$moment"
    cd "ZonalMD-$moment"
else 
    cd "$1"
    git pull
fi
mkdir -p build && cd build
cmake .. -DMD_FLEXIBLE_USE_MPI=ON -DCMAKE_BUILD_TYPE=Release -DMD_FLEXIBLE_FUNCTOR_AT_AUTOVEC=ON
make md-flexible -j 5
