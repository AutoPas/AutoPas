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
cd $HOME/code-server/scratch_space/tmp/transport

ls

tar -cvf all_code.tar *

bzip2 -9 all_code.tar 

cp all_code.tar.bz2 ../zipped

rm *.so

rm all_code.tar.bz2

## where do you want to send the files?

### THIS IS SUBJECT TO CHANGE FOR EXPERIMENTS ON DIFFERENT PLATFORMS
cd $HOME/code-server/scratch_space/tmp/zipped/

scp all_code.tar.bz2 tiwari@brood00:/hivehomes/tiwari/irs.1.0/decks/code/

ssh brood00 "tar -xjf /hivehomes/tiwari/irs.1.0/decks/code/all_code.tar.bz2 -C /hivehomes/tiwari/irs.1.0/decks/code/"
 
ssh brood00 "rm /hivehomes/tiwari/irs.1.0/decks/code/all_code.tar.bz2"
