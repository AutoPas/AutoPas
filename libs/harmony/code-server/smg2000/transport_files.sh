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
cd $HOME/scratch/new_code

tar -cvf all_code.tar *

bzip2 all_code.tar 

scp all_code.tar.bz2 tiwari@brood00:~/scratch/code/
rm *.so
rm all_code.tar.bz2

ssh tiwari@brood00 tar -xjf /hivehomes/tiwari/scratch/code/all_code.tar.bz2 -C /hivehomes/tiwari/scratch/code/

ssh tiwari@brood00 rm /hivehomes/tiwari/scratch/code/all_code.tar.bz2
