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
        iteration = 0

        regexConfigurationPrediction = '.* Traversal time prediction for +({.*}).*: *([0-9]+)'
        regexNoPrediction = '.* No traversal time prediction for +({.*})'
        regexIter = '.*Iteration +([0-9]+)'

        for line in file.readlines():
            if (match := re.search(regexIter, line)) is not None:
                iteration = match.group(1)
            elif (match := re.search(regexNoPrediction, line)) is not None:
                continue
            elif (match := re.search(regexConfigurationPrediction, line)) is not None:
                m = match.group(1)
                if match.group(1) in configurationPrediction:
                    configurationPrediction[match.group(1)].append((iteration, int(match.group(2))))
                else:
                    configurationPrediction[match.group(1)] = [(iteration, int(match.group(2)))]

    # create figure and define layout
    fig = go.Figure(
        layout=dict(
            showlegend=True,
            title_text=datafile,
            xaxis_title_text="Iteration",
            yaxis_title_text="Predicted time per Iteration",
        ),
    )

    # plotting predictions
    # probably a limited amount of configurations should be in one plot - probably like 5 or 10
    # maybe as user input if is should be in one or not
    # md-flexible test to test this script is not created yet so this is a prototype
    configsInPlot = 0
    for configuration in configurationPrediction:
        allPrediction = []
        allIteration = []

        for iteration, prediction in configurationPrediction[configuration]:
            allIteration.append(iteration)
            allPrediction.append(prediction)

        fig.add_trace(go.Scatter(x=allIteration, y=allPrediction, mode='lines+markers', name=configuration))

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