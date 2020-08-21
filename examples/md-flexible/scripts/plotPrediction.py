#!/usr/bin/python3

import sys
import plotly.graph_objects as go
import re
import os

# THIS SCRIPT NEEDS AT LEAST PYTHON 3.8.
# However, lesser version will probably fail due to invalid syntax instead of this assertion

# ---------------------------------------------- Input ----------------------------------------------
outputfileType = "none"
for arg in sys.argv[1:]:
    if "--help" in arg:
        print("Usage: ./plotPrediction.py OPTION FLAG path/To/mdFlex/std.out ...")
        print("Output options:\n "
              " prediction  - Shows the predictions for every configuration\n "
              " test        - Shows the tests for every configuration\n "
              "Flags:\n "
              " --png --jpeg --pdf - Plot is generated in the given file format")
        print("If no input is given the script does not work.")
        print("If no input is given the script does not work.")
        exit(0)
    elif "--png" in arg:
        outputfileType = "png"
    elif "--jpeg" in arg:
        outputfileType = "jpeg"
    elif "--pdf" in arg:
        outputfileType = "pdf"

# take all input files as source for a plot
if outputfileType != "none" and len(sys.argv) > 3:
    option = sys.argv[1]
    if option != "prediction" and option != "test" and sys.argv[2] != "--" + outputfileType:
        print("Error: Wrong input given! ./plotDiffPredictionTest.py --help to see what is needed.")
        sys.exit(-1)
    datafiles = sys.argv[3:]
elif len(sys.argv) > 2 and outputfileType == "none":
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

    # test if the file contains information that can be plotted
    # if len(configurationPrediction) == 0:
    #     print(datafile + ": No information could be extracted from this file!")
    #     continue

    # create figure and define layout
    text = "Test"
    if option == "prediction":
        text = "Predict"
    fig = go.Figure(
        layout=dict(
            showlegend=True,
            title_text=datafile,
            xaxis_title_text="Iteration",
            yaxis_title_text=text + "ed time per iteration",
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
    configurations = configurationPrediction
    if "test" == option:
        configurations = configurationTest
    i = 1
    for configuration in configurations:
        regexContainer = '.*Container: +(.*) , Cell.*'
        regexTraversal = '.*Traversal: +(.*) , Load.*'
        regexDataLayout = '.*Data Layout: +(.*) , Newton.*'
        regexNewton3 = '.*Newton 3: +(.*)}.*'

        container = ""
        if (match := re.search(regexContainer, configuration)) is not None:
            container = match.group(1)

        if "prediction" == option:
            allPrediction = []
            allIteration = []

            for iteration, prediction in configurationPrediction[configuration]:
                allIteration.append(iteration)
                allPrediction.append(prediction)

            if not shown[container]:
                fig.add_trace(
                    go.Scatter(x=allIteration, y=allPrediction, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]),
                               name=container))
                shown[container] = True
            else:
                fig.add_trace(
                    go.Scatter(x=allIteration, y=allPrediction, mode='lines+markers', legendgroup=container,
                               line=dict(color=colors[container]), name=container, showlegend=False))

        if "test" == option:
            allTest = []
            allIteration = []
            for iteration, test in configurationTest[configuration]:
                allIteration.append(iteration)
                allTest.append(test)

            if not shown[container]:
                fig.add_trace(go.Scatter(x=allIteration, y=allTest, mode='markers', legendgroup=container,
                                         line=dict(color=colors[container]),
                                         name=container))
                shown[container] = True
            else:
                fig.add_trace(go.Scatter(x=allIteration, y=allTest, mode='markers', legendgroup=container,
                                         line=dict(color=colors[container]),
                                         name=container, showlegend=False))

        i = i + 1

    if outputfileType == "none":
        fig.show()
    else:
        if not os.path.exists("images"):
            os.mkdir("images")

        fig.write_image("images/" + datafile + "." + outputfileType, scale=1.5)
