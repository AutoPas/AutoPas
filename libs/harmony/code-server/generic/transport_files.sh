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

appname=$1

cd $HOME/scratch/new_code_${appname}

tar -cvf all_code_${appname}.tar *

bzip2 all_code_${appname}.tar 

echo "done zipping"

# brood
scp all_code_${appname}.tar.bz2 rahulp@brood00:~/scratch/code/

########################## TRANSPORT ####################################


#hopper
rm -rf *.o *.c *.so
rm -rf all_code_${appname}.tar.bz2


# brood
ssh rahulp@brood00 tar -xjf /hivehomes/rahulp/scratch/code/all_code_${appname}.tar.bz2 -C /hivehomes/rahulp/scratch/code/
ssh rahulp@brood00 rm /hivehomes/rahulp/scratch/code/all_code_${appname}.tar.bz2



