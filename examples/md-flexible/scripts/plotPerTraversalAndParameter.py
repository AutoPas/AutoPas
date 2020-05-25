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
        print("Usage: ./plotTuning.py parameter=[number, size, density] [path/To//Output/*.out ...]. Meaning number = number of particles, size = boxSize and density = particle density")
        exit(0)
    elif parameterIdentifierString in arg:
        parameterArg = arg.split('=', 1)[1]
        containsParameterArg = True
    elif not containsParameterArg and not parameterIdentifierString in arg:
        print("Please specify the parameter for the x-axis: parameter=[number, size, density]")
        exit(0)

# take all input files as source for a plot
if len(sys.argv) > 1:
    if ".out" in arg:
        datafiles = sys.argv[2:]
    else:
        datadirs = sys.argv[2:]
        datadirs = list(datadirs)
        datadirs.sort(reverse=True)
        datafiles = os.listdir(datadirs[0])
        datafiles = filter(lambda s: s.endswith('.out'), datafiles)
        datafiles = map(lambda s: datadirs[0] + "/" + s, datafiles)
        datafiles = list(datafiles)

#    TODO: if nothing is given search for the latest test dir

# ---------------------------------------------- Script ---------------------------------------------

# values for plotting in the end
allTraversals = []
numberOfParticlesPerTraversal = []
boxSizePerTraversal = []
densityPerTraversal = []
timePerTravsersal = []
number = False
size = False
density = False
xAxisTitle = "empty"
xAxis = []

if parameterArg == 'number':
    number = True
    xAxisTitle = 'Number of particles'
elif parameterArg == 'size':
    size = True
    xAxisTitle = 'Boxsize'
elif parameterArg == 'density':
    density = True
    xAxisTitle = 'Density of particles'


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
    regexBoxMin = '.*box-min*'
    regexBoxMax = '.*box-max*'
    regexDensity = '.*particle-spacing*'
    regexTraversal = '.*traversal*'

    # parse file
    with open(datafile) as f:
        foundTraversal = False
        traversal = "noTraversal"

        counter = 0
        for line in f.readlines():
            if (match := re.search(regexTraversal, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                traversal = currentLine[0]
                if (not allTraversals.__contains__(traversal)):
                    allTraversals.append(traversal)
                    numberOfParticlesPerTraversal.append([])
                    boxSizePerTraversal.append([])
                    densityPerTraversal.append([])
                    timePerTravsersal.append([])
                foundTraversal = True
            elif number and (match := re.search(regexNumOfParticles, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')  # split content inside brackets and show as array
                numberOfParticles = numpy.prod(list(map(int, arrayOfCurrentLine)))  # calculate overall number of particles
                numberOfParticlesPerTraversal[allTraversals.index(traversal)].append(numberOfParticles)
            elif size and (match := re.search(regexBoxMax, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMax = list(map(float, arrayOfCurrentLine))
            elif size and (match := re.search(regexBoxMin, line)) is not None:
                currentLine = re.findall(r'\[(.*?)\]', line)  # get content inside the brackets
                arrayOfCurrentLine = currentLine[0].split(',')
                boxSizeListMin = list(map(float, arrayOfCurrentLine))
            elif density and (match := re.search(regexDensity, line)) is not None:
                currentLine = line.split(':', 1)[1]
                currentLine.strip()
                densityPerTraversal[allTraversals.index(traversal)].append(float(currentLine))
            elif (match := re.search(regexIterTook, line)) is not None:
                values.append(int(match.group(1)))  # append time needed to perform iteration
        foundTraversal = False

    timePerTravsersal[allTraversals.index(traversal)].append((sum(values)/len(values)))  # append mean time of this traversal
    if size:
        boxSize = numpy.prod([a - b for a, b in zip(boxSizeListMax, boxSizeListMin)])  # get size per dimension and then multiply to get overall boxsize
        boxSizePerTraversal[allTraversals.index(traversal)].append(boxSize)


# ---------------------------------------------- Plot ---------------------------------------------

# build the full rgb color space
numcolors=len(allTraversals)
distBetweenColors=int((6*256)/numcolors)
colorrange=[]
for g in range(0,256):
    colorrange.append('rgb(255, ' + str(g) + ',   0)')
for r in reversed(range(0,256)):
    colorrange.append('rgb(' + str(r) + ', 255,   0)')
for b in range(0,256):
    colorrange.append('rgb(  0, 255, ' + str(b) + ')')
for g in reversed(range(0,256)):
    colorrange.append('rgb(  0, ' + str(g) + ', 255)')
for r in range(0,256):
    colorrange.append('rgb(' + str(r) + ',   0, 255)')
for b in reversed(range(0,256)):
    colorrange.append('rgb(255,   0, ' + str(b) + ')')
# select equidistant colors
colorrange=colorrange[0::distBetweenColors]


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
    xAxis.sort()  # TODO: sort according to related time
    timePerTravsersal[allTraversals.index(t)].sort()
    fig.add_trace(go.Scatter(
        name=t,
        x=xAxis,
        y=timePerTravsersal[allTraversals.index(t)],
        mode="lines+markers",
        hovertext=t,
        marker=dict(
            color=allTraversals.index(t),
        ),
        line=dict(
            color=colorrange[allTraversals.index(t)],
        ),
        )
    )
fig.show()

