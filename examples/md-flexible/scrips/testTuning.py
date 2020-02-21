#!/bin/python3

import os
import re
import sys
import subprocess
from datetime import datetime

# ---------------------------------------- Global parameters and input ----------------------------------------

GREEN='\033[32m'
YELLOW='\033[93m'
RED='\033[31m'
ENDCOLOR='\033[0m'

# string to identify the argument identifying...
# the simulation binary
simulationIdentifierString='sim='
# the tuning strategy
tuningIdentifierString='tune='

# path to the testing binary
simulation='../../../cmake-build-manually/examples/md-flexible/md-flexible'
# placeholder for a globally set tuning strategy
tuningArg=[]
# list of input files or directories
configsDirs=[]

# parse special args
for arg in sys.argv[1:]:
    if "help" in arg:
        print("Usage: ./testTuning.py [sim=path/to/simulation/binary] [tune=tuningStrategy] [paths/to/yaml/files or/to/directories]")
        exit(0)
    elif simulationIdentifierString in arg:
        simulation=arg.split('=', 1)[1]
    elif tuningIdentifierString in arg:
        tuningArg=["--tuning-strategy", arg.split('=',1)[1]]
    else:
        # everything else is considered a path to inputs
        configsDirs.append(arg)

# sanitize parsed stuff
if not os.path.isfile(simulation):
    print("Not a file: " + simulation)
    exit(1)

simulation=os.path.abspath(simulation)

# default search directory for inputs:
if not configsDirs:
    configsDirs=['../input/testTuning/']

# directory for simulation output
outputDir="testTuning_"+datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

# ------------------------------------------------- Functions -------------------------------------------------

def getStringFromFile(filename, regex):
    with open(filename) as f:
        for line in f.readlines():
            # print(line, end='')
            match=re.search(regex, line)
            if match:
                return match.group(1)
                
def testScenario(yamlFile):
    scenarioName=yamlFile.split('/')[-1].split('.', 1)[0]
    print("Testing Scenario " + scenarioName)
    yamlFile=os.path.abspath(yamlFile)
    expected=getStringFromFile(yamlFile, '.* Expected Configuration +: +({.*})')
    command=[simulation, "--log-level", "debug" , "--no-end-config" , "--yaml-filename" , yamlFile] + tuningArg
    print(" ".join(command))
    outputFile=outputDir + '/' + scenarioName + '.out'
    with open(outputFile, 'w+') as outputLocation:
        subprocess.call(command, stdout=outputLocation)
    expectedTime=int(getStringFromFile(outputFile, expected + ".* Reduced value: ([0-9]+)"))
    selected=getStringFromFile(outputFile, '.* Selected Configuration +({.*})')
    selectedTime=int(getStringFromFile(outputFile, selected + ".* Reduced value: ([0-9]+)"))
    if expected == selected:
        print(GREEN + "Expected Configuration selected!" + ENDCOLOR)
    elif expectedTime > selectedTime:
        print(YELLOW + "Selected configuration faster than expected! Estimates maybe invalid for this hardware!" + ENDCOLOR)
    elif (expectedTime - selectedTime) < (expectedTime / 100):
        print(YELLOW + "Selected configuration less than 1% slower than expected one!" + ENDCOLOR)
    else:
        print(RED + "Inefficient configuration selected!" + ENDCOLOR)
    print("Selected: " + selected + " : " + str(selectedTime))
    print("Expected: " + expected + " : " + str(expectedTime))
    print()


# --------------------------------------------------- Script --------------------------------------------------

# create output directory
try:
    os.mkdir(outputDir)
except OSError:
    print("Could not create the output directory: " + outputDir)
    exit(2)

# iterate over all inputs and check for yaml files
for arg in configsDirs:
    if simulationIdentifierString in arg or tuningIdentifierString in arg:
        continue
    if os.path.isdir(arg):
        inputsFound=False
        for f in os.listdir(arg):
            if f.endswith('.yaml'):
                inputsFound=True
                testScenario(arg + f)
        if not inputsFound:
            print("No yaml files found in " + arg)
    elif os.path.isfile(arg):
        if arg.endswith('.yaml'):
            testScenario(arg)
        else:
            print("Not a yaml file: " + arg)
    else:
        print("Neither file nor folder: " + arg)

