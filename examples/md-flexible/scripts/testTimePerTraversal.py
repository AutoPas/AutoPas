#!/usr/bin/python3

import os
import re
import sys
import subprocess
from datetime import datetime

# ------------------------------------------------- Preprocessing Functions ------------------------------------

# extract all possible traversals from AllOptions
def getAllTraversals():
    regexTraversal = '.*traversal*'
    with open("AllOptions.yaml") as f:
        for line in f.readlines():
            if (match := re.search(regexTraversal, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                traversalList = currentLine[0].split(", ")
                for c in traversalList:
                    allTraversals.append(c)
                print(allTraversals)


# ---------------------------------------- Global parameters and input ----------------------------------------

# some definitions for colored text
GREEN='\033[32m'
YELLOW='\033[93m'
RED='\033[31m'
ENDCOLOR='\033[0m'

# Script variables can be overwritten:
# String to identify the argument identifying...
# the chosen traversal
traversalIdentifierString='traversal='

# the chosen amount of tests per traversal
testsIdentifierString='tests='

# path to the testing binary. This path is correct after cmake copied this script to the build dir.
simulation='./md-flexible'
# placeholder for used traversal
traversalArg=[]
containsTraversalArg = False
# placeholder for amount of tests
tests = 0
containsNumberOfTests = False
# list of input files or directories
configsDirs=[]
# list of all possible traversals
allTraversals=[]

# parse special args
for arg in sys.argv[1:]:
    if "help" in arg:
        print("Usage: ./testTimePerTraversal.py [traversal=chosenTraversal] [tests=iterations of tests] [paths/to/yaml/files or/to/directories]")
        print("Please use at least python 3.8.")
        exit(0)
    elif traversalIdentifierString in arg:
        traversalArg = ["--traversal", arg.split('=', 1)[1]]
        containsTraversalArg = True
    elif testsIdentifierString in arg:
        tests = int(arg.split('=', 1)[1])
        containsNumberOfTests = True
    else:
        # everything else is considered a path to inputs
        configsDirs.append(arg)
        if not containsTraversalArg:
            getAllTraversals()
        if not containsNumberOfTests:
            tests = 1

# sanitize parsed stuff
if not os.path.isfile(simulation):
    print("Not a file: " + simulation)
    exit(1)

simulation=os.path.abspath(simulation)

# default search directory for inputs:
if not configsDirs:
    configsDirs = ['../input/testTimePerBoxSizeOriginalSpacing/']

outputDir = ""
# directory for simulation
for arg in configsDirs:
    dirNames = arg.split('/')
    if os.path.isdir(arg):
        if dirNames[len(dirNames) - 1].__eq__(''):
            outputDir = outputDir + dirNames[len(dirNames)-2] + '_'
        else:
            outputDir = outputDir + dirNames[len(dirNames)-1] + '_'
outputDir = outputDir + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
outputDirName = "outputDir: " + outputDir
print(outputDirName)

# ------------------------------------------------- Functions -------------------------------------------------

# extracts first match of substring from a file
def getStringFromFile(filename, regex):
    with open(filename) as f:
        for line in f.readlines():
            match = re.search(regex, line)
            if match:
                return match.group(1)

# runs a given scenario and checks if the tuning result matches the expectation
def testScenario(yamlFile):
    # pretty print scenario name
    scenarioName=os.path.splitext(os.path.basename(yamlFile))[0]
    print("Testing Scenario " + scenarioName)
    yamlFile=os.path.abspath(yamlFile)

    # parse what configuration we expect to be optimal
    expected=getStringFromFile(yamlFile, '.* [eE]xpect.*({.*})')

    # build and execute command for simulation
    # append tuningArg list (if nothing is set this is empty)
    if len(allTraversals) > 0:
        for t in allTraversals:
            for test in range(0, tests):
                localTraversalArg = ["--traversal", t]
                command=[simulation, "--log-level", "debug", "--no-end-config", "--yaml-filename", yamlFile] + localTraversalArg
                print(" ".join(command))
                outputFile=os.path.join(outputDir, scenarioName + t + '_test' + str(test) + '.out')
                with open(outputFile, 'w+') as outputLocation:
                    subprocess.call(command, stdout=outputLocation, shell=False)

                # parse time of selected config
                selected=getStringFromFile(outputFile, '.* Selected Configuration +({.*})')
                selectedTime=int(getStringFromFile(outputFile, selected + ".* Reduced value: ([0-9]+)"))

                # print result
                print("Selected: " + selected + " : " + str(selectedTime))
                print()

    else:
        for test in range(0, tests):
            command=[simulation, "--log-level", "debug", "--no-end-config", "--yaml-filename", yamlFile] + traversalArg
            print(" ".join(command))
            outputFile=os.path.join(outputDir, scenarioName + traversalArg[1] + '_test' + str(test) + '.out')
            with open(outputFile, 'w+') as outputLocation:
                subprocess.call(command, stdout=outputLocation, shell=False)

            # parse time of selected config
            selected=getStringFromFile(outputFile, '.* Selected Configuration +({.*})')
            selectedTime=int(getStringFromFile(outputFile, selected + ".* Reduced value: ([0-9]+)"))

            # print result
            print("Selected: " + selected + " : " + str(selectedTime))
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
    if os.path.isdir(arg):
        inputsFound=False
        for f in os.listdir(arg):
            if f.endswith('.yaml'):
                inputsFound=True
                testScenario(os.path.join(arg, f))
        if not inputsFound:
            print("No yaml files found in " + arg)
    elif os.path.isfile(arg):
        if arg.endswith('.yaml'):
            testScenario(arg)
        else:
            print("Not a yaml file: " + arg)
    else:
        print("Neither file nor folder: " + arg)

