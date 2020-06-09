#!/usr/bin/python3

import sys
import re
from numpy import double

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion

# ---------------------------------------------- Input ----------------------------------------------

# help message
for arg in sys.argv[1:]:
    if "--help" in arg:
        print("Usage: ./compareOptimumConfiguration.py OPTION /to/mdFlexOutput1.out path/to/mdFlexOutput2.out "
              "[path/to/mdFlexOutputN.out]")
        print("Compare multiple md-flexible outputs for similarity between tuning results.")
        print(
            "The script need as input at least two files of output of md-flexible at log-level debug, where the first "
            "file is the base to which the other files are going to be compared to.")
        print("Optionally a float can be given and then the script will compare the evidence between tuning results.")
        print("If no or not enough input files are given the script does not work.")
        exit(0)

# take all input files as source for a plot
if len(sys.argv) > 2:
    try:
            relativeDifference = float(sys.argv[1])
            if len(sys.argv) > 3:
                datafiles = sys.argv[2:]
            else:
                print("Error: No or wrong input given! ./compareOptimumConfiguration.py --help to see what is needed.")
                sys.exit(-1)
    except ValueError:
        relativeDifference = 0
        datafiles = sys.argv[1:]
else:
    print("Error: No or wrong input given! ./compareOptimumConfiguration.py --help to see what is needed.")
    sys.exit(-1)

# ---------------------------------------------- Script ---------------------------------------------

comparisonBase = []
comparisonFiles = []

# collect data from files
for datafile in datafiles:

    # parse file
    with open(datafile, 'r') as file:

        # variables for data collection
        configurationsSelected = []
        configurationTest = {}

        regexSelectedConfiguration = '.* Selected Configuration +({.*})'
        regexCollectedTimes = '.* Collected times for +({.*})..*\[(.*)\].*: *([0-9]+)'

        for line in file.readlines():

            if (match := re.search(regexCollectedTimes, line)) is not None:
                configurationTest[match.group(1)] = int(match.group(3))
            elif (match := re.search(regexSelectedConfiguration, line)) is not None:
                configurationsSelected.append((match.group(1), configurationTest[match.group(1)]))

        if len(comparisonBase) == 0:
            comparisonBase = configurationsSelected
            # test if the base file contains information that can be compared
            if len(comparisonBase) == 0:
                print(datafiles[0] + ": Base file does not contain any information that can be compared!")
                exit(0)
        else:
            comparisonFiles.append(configurationsSelected)

# compare the files with the base file
j = 1
for fileToCompare in comparisonFiles:
    # test if the file contains information that can be compared
    if len(fileToCompare) == 0:
        print(datafiles[0] + " vs " + datafiles[j] + ":\n" + datafiles[j] + 
              "does not contain any information that can be compared!")
        continue

    loopInt = 0
    if len(comparisonBase) <= len(fileToCompare):
        loopInt = len(comparisonBase)
    else:
        loopInt = len(fileToCompare)
    countSameConfigSelected = 0
    countSameTimeSelected = 0

    i = 0
    while i < loopInt:
        [configBase, timeBase] = comparisonBase[i]
        [configToCompare, timeToCompare] = fileToCompare[i]
        if configBase == configToCompare:
            countSameConfigSelected += 1
        if timeToCompare < timeBase * relativeDifference + timeBase:
            countSameTimeSelected += 1

        i += 1

    print(datafiles[0] + " vs " + datafiles[j] + ":")
    print("Selected same configuration:                                  " + str(countSameConfigSelected) +
          "/" + str(loopInt) + " = " + str(double(countSameConfigSelected / loopInt * 100)) + "%")
    if relativeDifference != 0:
        print("Selected configuration that is within " + str(relativeDifference * 100) + "% of the reference: " + str(
            countSameTimeSelected) + "/" + str(loopInt) + " = " + str(double(countSameTimeSelected / loopInt * 100)) +
            "%")
    j += 1
