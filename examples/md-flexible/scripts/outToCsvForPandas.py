#!/usr/bin/python3

import os
import sys
import re
import numpy
import csv

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion
assert sys.version_info >= (3, 8)

# ---------------------------------------------- Input ----------------------------------------------

# Will be the name of the directory where the input files are stored
# if only one file is used, the filename will be the scenario name
nameOfScenario = None
setNameOfScenario = False

hardwareIdentifierString = 'hardware='
hardware = 'undefined'
lengthOfArguments = 1

# help message
for arg in sys.argv[1:]:
    if "--help" in arg:
        print(
            "Usage: ./outToCsvForPandas.py [hardware={name of hardware}] [path/To//Output/*.out ...]."
            "hardware: hardware the scenario was tested on")
        print("Please use at least python 3.8.")
        exit(0)
    elif hardwareIdentifierString in arg:
        hardware = arg.split('=', 1)[1]
        lengthOfArguments = 2

# take all input files as source for a plot
if len(sys.argv) > lengthOfArguments:
    if os.path.isdir(arg):
        datadirs = sys.argv[lengthOfArguments:]
        datadirs = list(datadirs)
        datadirs.sort(reverse=True)
        datafiles = os.listdir(datadirs[0])
        datafiles = filter(lambda s: s.endswith('.out'), datafiles)
        datafiles = map(lambda s: datadirs[0] + '/' + s, datafiles)
        datafiles = list(datafiles)
        nameOfDirs = datadirs[0].split('/')
        nameOfScenario = nameOfDirs[len(nameOfDirs) - 2]
    else:
        datafiles = sys.argv[lengthOfArguments:]
        setNameOfScenario = True

# -------------------------------------------- Functions --------------------------------------------

def parseConfigToDict(confStr):
    return {key.strip(): val.strip() for key, val in [pair.split(":") for pair in confStr[1:-1].split(',')]}


def appendForm(outerForm, appender):
    if outerForm is not None:
        innerForm = outerForm + '-' + appender
    else:
        innerForm = appender
    return innerForm


# ---------------------------------------------- Script ---------------------------------------------

# values for plotting in the end
configs = []
valuesErrorPlus = []
valuesErrorMinus = []
rowList = []

# gather data per file
for datafile in datafiles:
    # lists to gather data series
    values = []
    thisConfig = None
    finishedHeader = False

    regexConf = '.*Iterating with configuration: +({.*})'
    regexCellSizeFactor = '.* CellSizeFactor +({.*})'
    regexIterTook = '.* IteratePairwise took +([0-9]+) .*'
    regexNumOfParticles = '.*particles-per-dimension*'
    regexNumOfParticlesAbsolute = '.*numberOfParticles*'
    regexBoxMin = '.*box-min*'
    regexBoxMax = '.*box-max*'
    regexTraversal = '.*traversal*'
    regexSimulationStarted = 'Starting simulation...*'
    regexCubeUniform = '.*CubeUniform:*'
    regexCubeGrid = '.*CubeGrid:*'
    regexSphere = '.*Sphere:*'
    regexCubeGauss = '.*CubeGauss:*'
    regexHomogeneity = '.*Homogeneity*'

    if setNameOfScenario:
        datafileNames = datafile.split('/')
        nameOfScenario = datafileNames[len(datafileNames) - 1]

    with open(datafile) as f:
        currentDensity = 0.0
        foundTraversal = False
        traversal = "noTraversal"
        boxSizeListMin = []
        boxSizeListMax = []
        numberOfParticles = 0
        form = None
        container = None
        traversal = None
        cellSizeFactor = None
        loadEstimator = None
        dataLayout = None
        newton3 = None
        homogeneity = None

        lines = f.readlines()
        for line in lines:
            # parse header and footer
            if (match := re.search(regexSimulationStarted, line)) is not None:
                finishedHeader = True
            elif (match := re.search(regexNumOfParticles, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')  # split content inside brackets and show as array
                numberOfParticles = numpy.prod(
                    list(map(int, arrayOfCurrentLine)))  # calculate overall number of particles
            elif (match := re.search(regexNumOfParticlesAbsolute, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                numberOfParticles = currentLine
            elif (match := re.search(regexBoxMax, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMax = list(map(float, arrayOfCurrentLine))
            elif (match := re.search(regexBoxMin, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMin = list(map(float, arrayOfCurrentLine))
            elif (match := re.search(regexCubeUniform, line)) is not None:
                form = appendForm(form, 'CubeUniform')
            elif (match := re.search(regexCubeGrid, line)) is not None:
                form = appendForm(form, 'CubeGrid')
            elif (match := re.search(regexSphere, line)) is not None:
                form = appendForm(form, 'Sphere')
            elif (match := re.search(regexCubeGauss, line)) is not None:
                form = appendForm(form, 'CubeGauss')
            elif (match := re.search(regexHomogeneity, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                homogeneity = currentLine
        # parse Iterations
        for line in lines:
            if (match := re.search(regexConf, line)) is not None:
                thisConfig = parseConfigToDict(match.group(1))
                container = thisConfig.get("Container")
                traversal = thisConfig.get("Traversal")
                cellSizeFactor = thisConfig.get("CellSizeFactor")
                loadEstimator = thisConfig.get("Load Estimator")
                dataLayout = thisConfig.get("Data Layout")
                newton3 = thisConfig.get("Newton 3")
            elif (match := re.search(regexIterTook, line)) is not None:
                time = int(match.group(1))
                x = boxSizeListMax[0] - boxSizeListMin[0]
                y = boxSizeListMax[1] - boxSizeListMin[1]
                z = boxSizeListMax[2] - boxSizeListMin[2]
                rowList.append([traversal, container, newton3, dataLayout, cellSizeFactor, loadEstimator, form,
                                nameOfScenario, hardware, x, y, z, numberOfParticles, time, homogeneity])

# ---------------------------------------------- Write CSV ---------------------------------------------

with open('dataForPlotting.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(
        ["Traversal", "Container", "Newton3", "DataLayout", "CellsizeFactor", "LoadEstimator",
         "Form", "NameOfScenario", "Hardware", "x", "y", "z", "NumberOfParticles", "IterationTime", "Homogeneity"])
    writer.writerows(rowList)
