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

###########################################################################
# all the files required for code generation (suif files, for example)
# should be provided. For now, we will just look for required_files.dat
# in the appname directory.

# This script expects the following parameters:
# 1) Application name
# 2) Scratch directory on slave host.
# 3) Local hostname.
# 3) List of slave identifiers (eg: host_#)

appname=$1
scratch_dir=$2
local_host=$3
shift 3

for node in `echo $* | sort`
do
    machine_name=`echo $node | cut -d'_' -f 1`
    #echo $machine_name

    # make the directories
    # copy the relevant files
    dir_string=""
    if [ "$machine_name" = "$local_host" ];
    then
        if [ "$machine_name" != "$last_machine_name" ];
        then
            dir_string=$scratch_dir/new_code_$appname
            echo "Clearing/creating $dir_string"
            rm -rf $dir_string; mkdir -p $dir_string

            dir_string=$scratch_dir/zipped_$appname
            echo "Clearing/creating $dir_string"
            rm -rf $dir_string; mkdir -p $dir_string
        fi

	dir_string=$scratch_dir/${node}_$appname
        echo "Clearing/creating $dir_string"
        rm -rf $dir_string; mkdir -p $dir_string

	cp $appname/chill_script.$appname.sh $dir_string
	while read line
	do
	    # Comments in the required_files.dat contain a '#'
	    if [[ $line != *#* ]]; then
		cp $appname/$line $dir_string
	    fi
	done < "$appname/required_files.dat"

    else
	#echo "copying remote version"
        if [ "$machine_name" != "$last_machine_name" ]; then
            dir_string="$scratch_dir/new_code_$appname"
            echo "Clearing/creating $machine_name:$dir_string"
            ssh $machine_name "rm -rf $dir_string; mkdir -p $dir_string"

            dir_string="$scratch_dir/zipped_$appname"
            echo "Clearing/creating $machine_name:$dir_string"
            ssh $machine_name "rm -rf $dir_string; mkdir -p $dir_string"
        fi

        dir_string="$scratch_dir/${node}_$appname"
        echo "Clearing/creating $machine_name:$dir_string"
	ssh $machine_name "rm -rf $dir_string; mkdir -p $dir_string"

	scp $appname/chill_script.$appname.sh $machine_name:$dir_string
	while read line
	do
	    # Comments in the required_files.dat contain a '#'
            echo $line | grep "#" > /dev/null
	    if [ "$?" = "1" ]; then
                scp $appname/$line $machine_name:$dir_string
	    fi
	done < "$appname/required_files.dat"
    fi

    last_machine_name=$machine_name
done
