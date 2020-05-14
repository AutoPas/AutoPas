#!/usr/bin/python3

import sys
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import re

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion
assert sys.version_info >= (3, 8)

# ---------------------------------------------- Input ----------------------------------------------

# take all input files as source for a plot
if len(sys.argv) > 1:
    datafiles = sys.argv[1:]
else:
    print("Error: No input file given!")
    sys.exit(-1)

# -------------------------------------------- Functions --------------------------------------------


# ---------------------------------------------- Script ---------------------------------------------

# one plot for each given file
for datafile in datafiles:

    # parse file
    with open(datafile, 'r') as file:

        # variables for data collection
        configurationPrediction = {}
        configurationTest = {}
        configurationDiffTestPrediction = {}
        iteration = 0
        iterationBeginTuning = 0
        tuning = True

        regexConfigurationPrediction = '.* Traversal time prediction for +({.*}).*: *([0-9]+)'
        regexNoPrediction = '.* No traversal time prediction for +({.*})'
        regexCollectedTimes = '.* Collected times for +({.*})..*\[(.*)\].*: *([0-9]+)'
        regexIter = '.*Iteration +([0-9]+)'
        regexTuning = '.*Tuning: +([a-z]+)'

        for line in file.readlines():

            if (match := re.search(regexIter, line)) is not None:
                iteration = match.group(1)
            elif (match := re.search(regexNoPrediction, line)) is not None:
                continue
            elif (match := re.search(regexConfigurationPrediction, line)) is not None:
                    configurationPrediction[match.group(1)] = (iteration, int(match.group(2)))
            elif (match := re.search(regexCollectedTimes, line)) is not None:
                    configurationTest[match.group(1)] = (iteration, int(match.group(3)))
            elif (match := re.search(regexTuning, line)) is not None:
                if match.group(1).lower() == 'true':
                    if not tuning:
                        iterationBeginTuning = iteration
                    tuning = True
                elif match.group(1).lower() == 'false':
                    if tuning:
                        for configuration in configurationTest:
                            if configuration in configurationPrediction:
                                test = configurationTest[configuration]
                                prediction = configurationPrediction[configuration]
                                if (test[0] > iterationBeginTuning) and (prediction[0] == iterationBeginTuning):
                                    if configuration in configurationDiffTestPrediction:
                                        configurationDiffTestPrediction[configuration].append(
                                            (iterationBeginTuning, prediction[1] - test[1]))
                                    else:
                                        configurationDiffTestPrediction[configuration] = \
                                            [(iterationBeginTuning, prediction[1] - test[1])]
                    tuning = False

    # create figure and define layout
    fig = go.Figure(
        layout=dict(
            showlegend=True,
            title_text=datafile,
            xaxis_title_text="Iteration",
            yaxis_title_text="Predicted time - Tested time per Iteration",
        ),
    )

    # plotting predictions
    # probably a limited amount of configurations should be in one plot - probably like 5 or 10
    # maybe as user input if is should be in one or not
    # md-flexible test to test this script is not created yet so this is a prototype
    configsInPlot = 0
    for configuration in configurationDiffTestPrediction:
        allDiff = []
        allIteration = []

        for iteration, diff in configurationDiffTestPrediction[configuration]:
            allIteration.append(iteration)
            allDiff.append(diff)

        fig.add_trace(go.Scatter(x=allIteration, y=allDiff, mode='lines+markers', name=configuration))

        configsInPlot = configsInPlot + 1

        # if configsInPlot == 5:
        #    fig.show()
        #    configsInPlot = 0
        #    fig = go.Figure(
        #        layout=dict(
        #            showlegend=True,
        #            title_text=datafile,
        #            xaxis_title_text="Iteration",
        #            yaxis_title_text="Predicted time per Iteration",
        #        ),
        #    )
    fig.show()
