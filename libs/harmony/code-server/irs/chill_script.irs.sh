#!/bin/bash
#
# Copyright 2003-2016 Jeffrey K. Hollingsworth
#
# This file is part of Active Harmony.
#
# Active Harmony is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Active Harmony is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Active Harmony.  If not, see <http://www.gnu.org/licenses/>.
#

# command line arguments
#  the first argument indicates whether or not to produce new code
#  the second argument indicates the configuration is a primary or 
#   secondary (from speculative code generation)
#  the remaining arguments are for code generation.

primary_or_secondary=$1
produce_code=$2
#app_name=$3

# if the transformation parameters need to be mapped to another domain, for example
#   if we are looking at tile sizes that are multiples of 4, then fo the mapping
#   here.

TI=$(($3*1))
TJ=$(($4*1))
TK=$(($5*1))
UI=$6
UJ=$7

# exports: indicate where the CHiLL components and libraries are

export OMEGA_P=$HOME/omega
export SUIFHOME=$HOME/suifhome

export PATH=$HOME/chill/bin:${SUIFHOME}/i386-linux/bin/:${PATH}
export LD_LIBRARY_PATH=${SUIFHOME}/i386-linux/solib:${LD_LIBRARY_PATH}

# Working directory
WORK_DIR=$HOME/code-server/scratch_space/tmp/$8
TRANSPORT_DIR=$HOME/code-server/scratch_space/tmp/transport/
# cd into the work_directory

cd $WORK_DIR

# set the compilers and optimization flags
useICC=0
useGCC=1
GCC_COMMAND="gcc -O2 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse "
ICC_COMMAND="icc -O3 -xN -unroll0 "

# debugging
echo "Currently generating code for configuration: TI: $TI, TJ: $TJ, TK: $TK, UI: $UI, UJ: $UJ"

# remove previous variant of the kernel and the result

rm -rf rmatmult3_at.so rmatmult3_at_modified.c OUT__1__6229__.c temp.script peri.result.overall peri.result.iters

# output filename
out_file=rmatmult3_${TI}_${TJ}_${TK}_${UI}_${UJ}.so

if [ $produce_code -eq 1 ]; then
    # generate a new CHiLL script using the transformation parameters
    echo "
    source: rmatmult3_at.spd
    procedure: 0
    loop: 0

    TI=$TI
    TJ=$TJ
    TK=$TK
    UI=$UI
    UJ=$UJ


    permute([1,2,3])
    print
    known(kmax>=3)
    known(imax>=3)
    known(jmax>=3)
    tile(0,3,TI)
    tile(0,3,TJ)
    tile(0,3,TK)
    print
    
    # with all three loops unrolled
    #unroll(0,4,UJ)
    #unroll(0,5,UJ)
    #unroll(0,6,UI)    

    # with only UJ and UI
    unroll(0,5,UJ)
    unroll(0,6,UI)
   
    print
    " > temp.script; 

    # Run chill and the post-processor
    #chill temp.script;
    #s2c rmatmult3_at.lxf > rmatmult3_at_modified.c;
    ./chill.v.0.1.8 temp.script;
    ./s2c.v.0.1.8 rmatmult3_at.lxf > rmatmult3_at_modified.c;

    #cat rmatmult3_at_temp.c | sed 's/__out_argv\[15l\]/(int (\*)\[15\]) __out_argv\[15l\]/g' >  rmatmult3_at_temp_2.c
    #cat rmatmult3_at_temp_2.c | sed 's/extern void/extern \"C\" void/g' > rmatmult3_at_.c

    #rm rmatmult3_at_temp*.c
    # generate the .so library and copy it to the necessary directory

    if [ $useICC -eq 1 ]; then
        #echo "using icc"
        $ICC_COMMAND -c -fpic rmatmult3_at_.c
        $ICC_COMMAND -shared -lc -o  rmatmult3_at_.so rmatmult3_at_.o
	mv $out_file ../new_code/
	if [ $primary_or_secondary -eq 1 ]; then
	    # this is primary conf, move this to the transport dir
	    cp ../new_code/${out_file} $TRANSPORT_DIR
	    #cp $WORK_DIR/new_code/${out_file} $TRANSPORT_DIR
	fi
        echo "done creating the shared library"
    fi
    if [ $useGCC -eq 1 ]; then
        echo "using gcc"
        $GCC_COMMAND -c -fpic rmatmult3_at_modified.c
        $GCC_COMMAND -shared -lc -o  $out_file  rmatmult3_at_modified.o
	mv $out_file ../new_code/
	if [ $primary_or_secondary -eq 1 ]; then
	    # this is primary conf, move this to the transport dir
	    cp ../new_code/${out_file} $TRANSPORT_DIR
	fi
        echo "$6 is done creating the shared library"
    fi

fi

if [ $produce_code -eq 0 ]; then
    # simply move the out_file to transport directory
    cd $WORK_DIR
    cp ../new_code/${out_file} $TRANSPORT_DIR
fi
