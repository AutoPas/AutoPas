#!/usr/bin/python3

import sys
import plotly.graph_objects as go
import re
import os
from numpy import double

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion

# ---------------------------------------------- Input ----------------------------------------------
flag = "none"
for arg in sys.argv[1:]:
    if "--help" in arg:
        print("Usage: ./plotDiffPredictionTest.py OPTION FLAG path/To/mdFlex/std.out ...")
        print("Output options:\n "
              " relative - Shows the relative difference between prediction and test\n "
              " total    - Shows the total difference between prediction and test\n"
              "Flags:\n "
              " --png --jpeg --pdf - Plot is generated in the given file format")
        print("If no input is given the script does not work.")
        exit(0)
    elif "--png" in arg:
        flag = "png"
    elif "--jpeg" in arg:
        flag = "jpeg"
    elif "--pdf" in arg:
        flag = "pdf"

# take all input files as source for a plot
if flag != "none" and len(sys.argv) > 3:
    option = sys.argv[1]
    if option != "relative" and option != "total" and sys.argv[2] != "--" + flag:
        print("Error: Wrong input given! ./plotDiffPredictionTest.py --help to see what is needed.")
        sys.exit(-1)
    datafiles = sys.argv[3:]
elif len(sys.argv) > 2 and flag == "none":
    option = sys.argv[1]
    if option != "relative" and option != "total":
        print("Error: Wrong input given! ./plotDiffPredictionTest.py --help to see what is needed.")
        sys.exit(-1)
    datafiles = sys.argv[2:]
else:
    print("Error: No input file given! ./plotDiffPredictionTest.py --help to see what is needed.")
    sys.exit(-1)

# ---------------------------------------------- Script ---------------------------------------------

# one plot for each given file
for datafile in datafiles:

    # parse file
    with open(datafile, 'r') as file:

        # variables for data collection
        configurationPrediction = {}
        configurationTest = {}
        configurationDiffTestPredictionTotal = {}
        configurationDiffTestPredictionRelative = {}
        iteration = 0
        iterationBeginTuning = 0
        tuning = True

        regexConfigurationPrediction = '.* Traversal time prediction for +({.*}).*: *([0-9]+)'
        regexNoPrediction = '.* No traversal time prediction for +({.*})'
        regexCollectedTimes = '.* Collected times for +({.*})..*\[(.*)\].*: *([0-9]+)'
        regexIter = '.*Iteration +([0-9]+)'
        regexTuning = '.*tuning: +([a-z]+)'

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
                                    if configuration in configurationDiffTestPredictionTotal:
                                        configurationDiffTestPredictionTotal[configuration].append(
                                            (iterationBeginTuning, prediction[1] - test[1]))
                                        if prediction[1] == 0:
                                            configurationDiffTestPredictionRelative[configuration].append(
                                                (iterationBeginTuning, double(1 / test[1])))
                                        else:
                                            configurationDiffTestPredictionRelative[configuration].append(
                                                (iterationBeginTuning, double(prediction[1] / test[1])))
                                    else:
                                        configurationDiffTestPredictionTotal[configuration] = \
                                            [(iterationBeginTuning, prediction[1] - test[1])]
                                        if prediction[1] == 0:
                                            configurationDiffTestPredictionRelative[configuration] = \
                                                [(iterationBeginTuning, double(1 / test[1]))]
                                        else:
                                            configurationDiffTestPredictionRelative[configuration] = \
                                                [(iterationBeginTuning, double(prediction[1] / test[1]))]
                    tuning = False

    # test if the file contains information that can be plotted
    if len(configurationDiffTestPredictionTotal) == 0 and ("total" == option or "both == option"):
        print(datafile + ": No information could be extracted from this file!")
        continue
    if len(configurationDiffTestPredictionRelative) == 0 and "relative" == option:
        print(datafile + ": No information could be extracted from this file!")
        continue

    # create figure and define layout
    fig = go.Figure(
        layout=dict(
            showlegend=True,
            title_text=datafile,
            xaxis_title_text="Iteration",
            yaxis_title_text=option+" difference between predicted time and tested time per iteration",
        ),
    )

    colors = {
        "DirectSum": '#808000',  # olive
        "LinkedCells": '#FF0000',  # red
        "VerletLists": '#008000',  # green
        "VerletListsCells": '#0000FF',  # blue
        "VerletClusterLists": '#4B0082',  # indigo
        "VarVerletListsAsBuild": '#FFA500',  # orange
        "VerletClusterCells": "#90EE90"  # lightgreen
    }

    shown = {
        "DirectSum": False,
        "LinkedCells": False,
        "VerletLists": False,
        "VerletListsCells": False,
        "VerletClusterLists": False,
        "VarVerletListsAsBuild": False,
        "VerletClusterCells": False
    }

    # plotting predictions
    i = 1
    for configuration in configurationDiffTestPredictionTotal:
        regexContainer = '.*Container: +(.*) , Cell.*'
        regexTraversal = '.*Traversal: +(.*) , Load.*'
        regexDataLayout = '.*Data Layout: +(.*) , Newton.*'
        regexNewton3 = '.*Newton 3: +(.*)}.*'

        container = ""
        if (match := re.search(regexContainer, configuration)) is not None:
            container = match.group(1)

        if "total" == option:
            allDiffTotal = []
            allIterationTotal = []
            for iteration, diff in configurationDiffTestPredictionTotal[configuration]:
                allIterationTotal.append(iteration)
                allDiffTotal.append(diff)

            if not shown[container]:
                fig.add_trace(
                    go.Scatter(x=allIterationTotal, y=allDiffTotal, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]), name=container))
                shown[container] = True
            else:
                fig.add_trace(
                    go.Scatter(x=allIterationTotal, y=allDiffTotal, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]), name=container, showlegend=False))

        if "relative" == option:
            allDiffRelative = []
            allIterationRelative = []
            for iteration, diff in configurationDiffTestPredictionRelative[configuration]:
                allIterationRelative.append(iteration)
                allDiffRelative.append(diff)

            if not shown[container]:
                fig.add_trace(
                    go.Scatter(x=allIterationRelative, y=allDiffRelative, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]), name=container))
                shown[container] = True
            else:
                fig.add_trace(
                    go.Scatter(x=allIterationRelative, y=allDiffRelative, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]), name=container, showlegend=False))
                fig.update_yaxes(range=[0, 2])
        i = i + 1

    if flag == "none":
        fig.show()
    else:
        if not os.path.exists("images"):
            os.mkdir("images")

        fig.write_image("images/" + datafile + "." + flag, scale=1.5)
