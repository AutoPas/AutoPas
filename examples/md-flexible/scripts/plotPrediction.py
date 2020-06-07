#!/usr/bin/python3

import sys
import plotly.graph_objects as go
import re

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion

# ---------------------------------------------- Input ----------------------------------------------
for arg in sys.argv[1:]:
    if "--help" in arg:
        print("Usage: ./plotPrediction.py OPTION path/To/mdFlex/std.out ...")
        print("Output options:\n "
              " prediction  - Shows the predictions for every configuration\n "
              " test        - Shows the predictions and tests for every configuration")
        print("If no input is given the script does not work.")
        exit(0)

# take all input files as source for a plot
if len(sys.argv) > 2:
    option = sys.argv[1]
    if option != "prediction" and option != "test":
        print("Error: Wrong input given! ./plotPrediction.py --help to see what is needed.")
        sys.exit(-1)
    datafiles = sys.argv[2:]
else:
    print("Error: No input file given! ./plotPrediction.py --help to see what is needed.")
    sys.exit(-1)

# ---------------------------------------------- Script ---------------------------------------------

# one plot for each given file
for datafile in datafiles:

    # parse file
    with open(datafile, 'r') as file:

        # variables for data collection
        configurationPrediction = {}
        configurationTest = {}
        iteration = 0

        regexConfigurationPrediction = '.* Traversal time prediction for +({.*}).*: *([0-9]+)'
        regexNoPrediction = '.* No traversal time prediction for +({.*})'
        regexCollectedTimes = '.* Collected times for +({.*})..*\[(.*)\].*: *([0-9]+)'
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
            elif (match := re.search(regexCollectedTimes, line)) is not None:
                if match.group(1) in configurationTest:
                    configurationTest[match.group(1)].append((iteration, int(match.group(3))))
                else:
                    configurationTest[match.group(1)] = [(iteration, int(match.group(3)))]

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
    for configuration in configurationPrediction:
        allPrediction = []
        allIteration = []

        for iteration, prediction in configurationPrediction[configuration]:
            allIteration.append(iteration)
            allPrediction.append(prediction)

        fig.add_trace(go.Scatter(x=allIteration, y=allPrediction, mode='lines+markers', name=configuration))

        if "test" == option:
            # do not know if this is the right solution testing needed and some research on other options.
            allTest = []
            allIteration = []
            for iteration, test in configurationTest[configuration]:
                allIteration.append(iteration)
                allPrediction.append(test)

            fig.add_trace(go.Scatter(x=allIteration, y=allPrediction, mode='markers', name=configuration))

    fig.show()
