#!/bin/sh
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
i=$1
j=$2
k=$3


echo "

    source: rmatmult3_at.spd
    procedure: 0
    loop: 0
    
    UI=$i
    UJ=$j
    UK=$k
    

    known(kmax>=3)
    known(imax>=3)
    known(jmax>=3)
    permute([1,2,3])

    tile(0,2,20)
    tile(0,2,20)
    tile(0,5,20)
    #print
    
    # with all three loops unrolled
    unroll(0,4,UK)
    unroll(0,5,UJ)
    unroll(0,6,UI)    

    # with only UJ and UI
    #unroll(0,5,UJ)
    #unroll(0,6,UI)

    

    print
" > temp.script;
./chill.v.0.1.8 temp.script;
./s2c.v.0.1.8 rmatmult3_at.lxf > rmatmult3_at_modified.c;
GCC_COMMAND="gcc -O2 -march=pentium4 -mmmx -msse -msse2 -mfpmath=sse ";


$GCC_COMMAND -c -fpic rmatmult3_at_modified.c;
#$GCC_COMMAND -shared -lc -o  $out_file  rmatmult3_at_modified.o;
du -b rmatmult3_at_modified.o > du_${i}_${j}_${k};

#out_file=rmatmult3_${i}_${j}_${k}.so;

#$GCC_COMMAND -shared -lc -o  $out_file  rmatmult3_at_modified.o;

