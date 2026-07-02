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

# generate a new CHiLL script using the transformation parameters

# this script gets the code transformation parameters as an argument.

# This generic (MM) example has five parameters: TI, TJ, TK, UI and UJ.

TI=$1
TJ=$2
TK=$3
UI=$4
UJ=$5

# name of the out_file
__out_file_prefix=${file_prefix}_${TI}_${TJ}_${TK}_${UI}_${UJ}

# message
echo "Currently generating code for configuration: TI: $TI, TJ: $TJ, TK: $TK, UI: $UI, UJ: $UJ"

# define and export the __out_file env variable. This is read by the top-level
#  chill_script_generic.(remote)+.sh

echo "
    source: ${file_prefix}.sp2
    format: suif
    procedure: 0
    loop: 0

    TI=$TI
    TJ=$TJ
    TK=$TK
    UI=$UI
    UJ=$UJ


    permute([3,1,2])

    tile(0,2,TJ)
    tile(0,2,TI)
    tile(0,5,TK)
    
    # unroll I and J
    unroll(0,4,UI)
    unroll(0,5,UJ)
   
    " > temp.script; 

