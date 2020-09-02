#!/usr/bin/python3

import os
import sys
import plotly.graph_objects as go
import re
import numpy

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion
assert sys.version_info >= (3, 8)

# ---------------------------------------------- Input ----------------------------------------------

parameterIdentifierString = 'parameter='
parameterArg = []
containsParameterArg = False

# help message
for arg in sys.argv[1:]:
    if "--help" in arg:
        print(
            "Usage: ./plotPerTraversalAndParameter.py parameter=[number, size, density, homogeneity] [path/To//Output/*.out ...]. Meaning number = number of particles, size = boxSize and density = particle density")
        print("Please use at least python 3.8.")
        exit(0)
    elif parameterIdentifierString in arg:
        parameterArg = arg.split('=', 1)[1]
        containsParameterArg = True
    elif not containsParameterArg and not parameterIdentifierString in arg:
        print("Please specify the parameter for the x-axis: parameter=[number, size, density]")
        exit(0)

# take all input files as source for a plot
if len(sys.argv) > 1:
    if os.path.isdir(arg):
        datadirs = sys.argv[2:]
        datadirs = list(datadirs)
        datadirs.sort(reverse=True)
        datafiles = os.listdir(datadirs[0])
        datafiles = filter(lambda s: s.endswith('.out'), datafiles)
        datafiles = map(lambda s: datadirs[0] + '/' + s, datafiles)
        datafiles = list(datafiles)
    else:
        datafiles = sys.argv[2:]

# ---------------------------------------------- Script ---------------------------------------------

# values for plotting in the end
allTraversals = []
allContainers = [[]] * 7
numberOfParticlesPerTraversal = []
boxSizePerTraversal = []
densityPerTraversal = []
homogeneityPerTraversal = []
timePerTravsersal = []
number = False
size = False
density = False
homogeneity = False
xAxisTitle = "empty"
xAxis = []
meanXaxis = []
meanYaxis = []
valuesErrorPlus = []
valuesErrorMinus = []

if parameterArg == 'number':
    number = True
    xAxisTitle = 'Number of particles'
elif parameterArg == 'size':
    size = True
    xAxisTitle = 'Boxsize'
elif parameterArg == 'density':
    density = True
    xAxisTitle = 'Density of particles'
elif parameterArg == 'homogeneity':
    homogeneity = True
    xAxisTitle = 'Standard Deviation of Homogeneity'

allContainers[0] = ["directSum"]  # directSum
allContainers[1] = ["c01", "c08", "c18", "sliced", "c01-combined-SoA", "c04", "c04SoA", "c04HCP"]  # linkedCells
allContainers[2] = ["verlet-clusters", "verlet-clusters-coloring", "verlet-clusters-static"]  # verletClusterLists
allContainers[3] = ["verlet-sliced", "verlet-c18", "verlet-c01"]  # verletListsCells
allContainers[4] = ["verlet-cluster-cells"]  # verletClusterCells
allContainers[5] = ["var-verlet-lists-as-build"]  # varVerletListsAsBuild
allContainers[6] = ["verlet-lists"]  # verletLists

# gather data per file
for datafile in datafiles:

    # lists to gather data series
    values = []
    numberOfParticles = []
    boxSize = []
    boxSizeListMin = []
    boxSizeListMax = []
    traversal = []

    regexIterTook = '.* IteratePairwise took +([0-9]+) .*'
    regexNumOfParticles = '.*particles-per-dimension*'
    regexNumOfParticlesAbsolute = '.*numberOfParticles*'
    regexBoxMin = '.*box-min*'
    regexBoxMax = '.*box-max*'
    regexTraversal = '.*traversal*'
    regexHomogeneity = '.*Homogeneity*'

    # parse file
    with open(datafile) as f:
        currentDensity = 0.0
        foundTraversal = False
        traversal = "noTraversal"

        counter = 0
        for line in f.readlines():
            if (match := re.search(regexTraversal, line)) is not None:
                currentLine = re.findall('\[(.*?)\]', line)  # get content inside the brackets
                if not foundTraversal:
                    traversal = currentLine[0]
                    if (not allTraversals.__contains__(traversal)):
                        allTraversals.append(traversal)
                        numberOfParticlesPerTraversal.append([])
                        boxSizePerTraversal.append([])
                        densityPerTraversal.append([])
                        timePerTravsersal.append([])
                        homogeneityPerTraversal.append([])
                    foundTraversal = True
            elif (number or density) and (match := re.search(regexNumOfParticles, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')  # split content inside brackets and show as array
                numberOfParticles = numpy.prod(
                    list(map(int, arrayOfCurrentLine)))  # calculate overall number of particles
                numberOfParticlesPerTraversal[allTraversals.index(traversal)].append(numberOfParticles)
                currentDensity += numberOfParticles
            elif number and (match := re.search(regexNumOfParticlesAbsolute, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                numberOfParticlesPerTraversal[allTraversals.index(traversal)].append(float(currentLine))
            elif (size or density) and (match := re.search(regexBoxMax, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMax = list(map(float, arrayOfCurrentLine))
            elif (size or density) and (match := re.search(regexBoxMin, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMin = list(map(float, arrayOfCurrentLine))
            elif density and (match := re.search(regexNumOfParticlesAbsolute, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                currentDensity = float(currentLine)
            elif homogeneity and (match := re.search(regexHomogeneity, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                homogeneityPerTraversal[allTraversals.index(traversal)].append(float(currentLine))
            elif (match := re.search(regexIterTook, line)) is not None:
                values.append(int(match.group(1)))  # append time needed to perform iteration
        foundTraversal = False

    timePerTravsersal[allTraversals.index(traversal)].append(
        (sum(values) / len(values)))  # append mean time of this traversal
    if size or density:
        boxSize = numpy.prod([a - b for a, b in zip(boxSizeListMax,
                                                    boxSizeListMin)])  # get size per dimension and then multiply to get overall boxsize
        boxSizePerTraversal[allTraversals.index(traversal)].append(boxSize)
        if density:
            currentDensity = currentDensity / boxSize
            densityPerTraversal[allTraversals.index(traversal)].append(currentDensity)


# ---------------------------------------------- Functions ----------------------------------------

def getContainerByTraversal(traversal):
    for c in allContainers:
        if c.__contains__(traversal):
            return c


def getTraversalIndexInContainer(traversal):
    for c in allContainers:
        if c.__contains__(traversal):
            return c.index(traversal)
    return 0


def calculateMeans(xAxis, yAxis):
    global meanYaxis
    global meanXaxis
    global valuesErrorMinus
    global valuesErrorPlus
    meanXaxis = []
    meanYaxis = []
    valuesErrorMinus = []
    valuesErrorPlus = []
    counter = 0
    for x in xAxis:
        if meanXaxis.__contains__(x):
            # newIndex means the place in meanYaxis where the values corresponding to the meanXaxis are currently stored
            newIndex = meanXaxis.index(x)
            # oldIndex references where the Y value was found within the yAxis array
            oldIndex = counter
            meanYaxis[newIndex] += yAxis[oldIndex]
            # update the ErrorMinus and Errorplus arrays if a smaller or larger value was found
            if valuesErrorMinus[newIndex] > yAxis[oldIndex]:
                valuesErrorMinus[newIndex] = yAxis[oldIndex]
            if valuesErrorPlus[newIndex] < yAxis[oldIndex]:
                valuesErrorPlus[newIndex] = yAxis[oldIndex]
            counter += 1
        else:
            # append x in x-axis and y ind y-axis as well as for the highest and lowest y value so far
            meanXaxis.append(x)
            meanYaxis.append(yAxis[xAxis.index(x)])
            valuesErrorMinus.append(yAxis[xAxis.index(x)])
            valuesErrorPlus.append(yAxis[xAxis.index(x)])
            counter += 1
    divisor = len(yAxis) / len(meanYaxis)
    meanYaxis = [y / divisor for y in meanYaxis]

    meanYaxis = [y for _, y in sorted(zip(meanXaxis, meanYaxis), key=lambda pair: pair[0])]
    valuesErrorPlus = [y for _, y in sorted(zip(meanXaxis, valuesErrorPlus), key=lambda pair: pair[0])]
    valuesErrorMinus = [y for _, y in sorted(zip(meanXaxis, valuesErrorMinus), key=lambda pair: pair[0])]
    meanXaxis.sort()


# ---------------------------------------------- Plot ---------------------------------------------

# build the full rgb color space
colorrange = ['rgb(  0, 0, 0)', 'rgb(0, 100,   100)', 'rgb(255, 219,   0)', 'rgb(255, 0,   0)', 'rgb(  0, 147, 255)',
              'rgb(71, 0, 255)', 'rgb(73, 255,   0)', 'rgb(255,   0, 221)']
symbolrange = ['circle', 'square', 'cross', 'star', 'diamond', 'star-diamond', 'bowtie', 'hourglass', 'triangle-up',
               'triangle-down', 'triangle-left', 'triangle-right']

# create figure and define layout
fig = go.Figure(
    layout=dict(
        title_text="Plot per " + xAxisTitle.lower() + " and per traversal",
        xaxis_title_text=xAxisTitle,
        yaxis_title_text="Time per Iteration [ns]",
    ),
)

# add data
for t in allTraversals:
    xAxis = []
    if number:
        xAxis = numberOfParticlesPerTraversal[allTraversals.index(t)]
    elif size:
        xAxis = boxSizePerTraversal[allTraversals.index(t)]
    elif density:
        xAxis = densityPerTraversal[allTraversals.index(t)]
    elif homogeneity:
        xAxis = homogeneityPerTraversal[allTraversals.index(t)]
    calculateMeans(xAxis, timePerTravsersal[allTraversals.index(t)])
    upperErrorBound = [x1 - x2 for (x1, x2) in zip(valuesErrorPlus, meanYaxis)]
    lowerErrorBound = [y1 - y2 for (y1, y2) in zip(meanYaxis, valuesErrorMinus)]
    fig.add_trace(go.Scatter(
        name=t,
        x=meanXaxis,
        y=meanYaxis,
        mode="lines+markers",
        error_y=dict(
            type='data',
            symmetric=False,
            array=upperErrorBound,
            arrayminus=lowerErrorBound,
        ),
        hovertext=t,
        marker=dict(
            color=colorrange[allContainers.index(getContainerByTraversal(t))],
            symbol=symbolrange[getTraversalIndexInContainer(t)],
            size=10,
        ),
        line=dict(
            color=colorrange[allContainers.index(getContainerByTraversal(t))]
        ),
    )
    )
fig.show()
