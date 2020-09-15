#!/usr/bin/python3

import os
import sys
import plotly.graph_objects as go
import re

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion
assert sys.version_info >= (3,8)

# ---------------------------------------------- Input ----------------------------------------------

# help message
for arg in sys.argv[1:]:
    if "--help" in arg:
        print("Usage: ./plotTuning.py [path/To/mdFlex/std.out ...] or [path/To/mdFlex/directoryWithOutput]")
        print("If no input is given the script looks for the latest testTuning directory in the current directory.")
        exit(0)

# take all input files as source for a plot
if len(sys.argv) > 1:
    if os.path.isdir(arg):
        datadirs = sys.argv[1:]
        datadirs = list(datadirs)
        datadirs.sort(reverse=True)
        datafiles = os.listdir(datadirs[0])
        datafiles = filter(lambda s: s.endswith('.out'), datafiles)
        datafiles = map(lambda s: datadirs[0] + '/' + s, datafiles)
        datafiles = list(datafiles)
    else:
        datafiles = sys.argv[1:]
else:
    # if nothing is given search for the latest test dir
    datadirs=os.listdir("./")
    datadirs=list(filter(lambda s : s.startswith('testTuning_'), datadirs))
    datadirs.sort(reverse = True)
    datafiles=os.listdir(datadirs[0])
    datafiles=filter(lambda s : s.endswith('.out'), datafiles)
    datafiles=map(lambda s : datadirs[0] + "/" + s, datafiles)
    datafiles=list(datafiles)

# -------------------------------------------- Functions --------------------------------------------

def parseConfigToDict(confStr):
    return {key.strip():val.strip() for key,val in [pair.split(":") for pair in confStr[1:-1].split(',')]}

# ---------------------------------------------- Script ---------------------------------------------

# one plot for each given file
for datafile in datafiles:

    # lists to gather data series
    iterationNr = []
    configs = []
    valuesErrorPlus = []
    valuesErrorMinus = []
    values = []

    regexCollectedTimes='.* Collected times for +({.*}).*\[(.*)\].*: *([0-9]+)'
    regexSelectedConf='.* Selected Configuration +({.*})'
    regexIterTook='.* IteratePairwise took +([0-9]+) .*'
    regexIter='.*Iteration +([0-9]+)'
    regexTuning='.*Tuning: +([a-z]+)'

    # parse file
    with open(datafile) as f:
        tuning=True
        selectedConf=None
        thisConfig=None
        thisIteration=0
        thisValue=0
        thisValuesErrorPlus=0
        thisValuesErrorMinus=0

        counter=0
        for line in f.readlines():
            if (match := re.search(regexIter, line)) is not None:
                thisIteration=match.group(1)
            elif (match := re.search(regexCollectedTimes, line)) is not None:
                thisTimes=[int(s) for s in match.group(2).split()]
                thisConfig=parseConfigToDict(match.group(1))
                thisValue=int(match.group(3))
                thisValuesErrorPlus=max(thisTimes)-thisValue
                thisValuesErrorMinus=thisValue-min(thisTimes)
            elif (match := re.search(regexSelectedConf, line)) is not None:
                selectedConf=parseConfigToDict(match.group(1))
                thisConfig=selectedConf
            elif not tuning and (match := re.search(regexIterTook, line)) is not None:
                thisConfig=selectedConf
                thisValue=match.group(1)
                thisValuesErrorPlus=0
                thisValuesErrorMinus=0
            elif (match := re.search(regexTuning, line)) is not None:
                if match.group(1).lower() == 'true':
                    tuning=True
                elif match.group(1).lower() == 'false':
                    tuning=False
                else:
                    print("Could not pase bool of tuning line!")
                    exit(1)
                # ASSUMPTION Tuning indicator is the last relevant output of an iteration
                if not tuning or thisConfig is not None:
                    iterationNr.append(thisIteration)
                    configs.append(thisConfig)
                    thisConfig=None
                    values.append(thisValue)
                    valuesErrorPlus.append(thisValuesErrorPlus)
                    valuesErrorMinus.append(thisValuesErrorMinus)


    # sanity checks:
    assert len(iterationNr) == len(configs)
    assert len(iterationNr) == len(values)
    assert len(iterationNr) == len(valuesErrorPlus)
    assert len(iterationNr) == len(valuesErrorMinus)

    # this set comprehension eliminates duplicates
    allContainers=list({c["Container"] for c in configs})
    allContainers.sort()

    # build the full rgb color space
    numcolors=len(allContainers)
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

    # build custom color scale by assigning each container a color
    colorscale=[]
    sectionWidth=1.0/len(allContainers)
    colorWidth=(256**3-1) / len(allContainers)
    for c in allContainers:
        containerIndex=allContainers.index(c)
        sectionStart=sectionWidth*containerIndex
        sectionStop=sectionWidth*(containerIndex+1)
        sectionColor=colorrange[allContainers.index(c)]
        colorscale.append([sectionStart,sectionColor])
        colorscale.append([sectionStop,sectionColor])

    # create figure and define layout
    fig = go.Figure(
            layout=dict(
                showlegend=False,
                title_text=datafile,
                xaxis_title_text="Iteration",
                yaxis_title_text="Time per Iteration [ns]",
                ),
            )

    # add data
    fig.add_trace(go.Scatter(
        x=iterationNr,
        y=values,
        mode="lines+markers",
        error_y=dict(
            type='data',
            symmetric=False,
            array=valuesErrorPlus,
            arrayminus=valuesErrorMinus,
            ),
        hovertext=configs,
        marker=dict(
            color=[allContainers.index(c["Container"]) for c in configs],
            showscale=True,
            colorscale=colorscale,
            colorbar=dict(
                tickmode="array",
                # center the labels of the colorscale
                tickvals=[(1-1/len(allContainers))*(x+.5) for x in range(0, len(allContainers))],
                ticktext=allContainers,
                ),
            ),
        line=dict(
            color='#AAAAAA',
            dash='dash',
            ),
        ))
    fig.show()
