# MD-Flex Scrips

This folder contains scripts for executing tests, measurements, and visualizations.

## Tuning tests

### testTuning.py

Requirements:
* python3 (tested with 3.6.9)

Tests if md-flexible selects the expected configuration for given scenarios. Examples can be found in:
AutoPas/examples/md-flexible/input/testTuning

### plotTuning.py

Requirements:
* python3.8 (tested with 3.8.2)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.5.1)

Creates an interactive plot for each md-flex output file given. If no input is given it tries to find the last output folder of testTuning.py and plot everything in there.

### plotPrediction.py

Requirements:
* python3 (tested with Python 3.8.3rc1)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.7.1)

Creates a plot of the predictions or the prediction and the tested times for each configuration for each md-flex output file given depending on the chosen option in the input.

### plotDiffPredictionTest.py

Requirements:
* python3 (tested with Python 3.8.3rc1)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.7.1)

Creates a plot of the total or the relative difference between the predictions and the tested times for each configuration for each md-flex output file given or both depending on the chosen option in the input.

### compareOptimumConfiguration.py

Requirements:
* python3 (tested with Python 3.8.3rc1)

Compares the tuning behavior of multiple runs of md-flexible. Specifically it compares the configuration which is selected as optimum after each tuning phase with a base file and every other file given in the input and prints the percentage of the alignment.
If the base file has fewer tuning phases than the compared file, than the rest of the compared file is not going to be compared.

### extractClusterGraph.py

Requirements:
* python3.8 (tested with 3.8.2)

Creates graphs for a md-flex output file of a bayesian-cluster run. Graphs consist of two csv files: one containing information for all nodes(clusters) and one for all edges(weights between clusters).
extractClusterGraphSeparate.py creates one file-pair for each graph found in the md-flex output. extractClusterGraphCombined.py combines all found graphs into one file-pair.
Additionaly a md-flex output of a fullSearch run can be provided. This appends the sampled runtimes to the corresponding nodes.

## Compare Traversals per time

Enables you to compare several traversals per time and number of particles or domain size or density.
In the plot you will find the time as the x-axis and the other parameter as y-axis.

Requirements:
* python3 (tested with 3.8.2)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.7.1)

### testTimePerTraversal.py
Test several traversals at once.

Examples can be found in AutoPas/examples/md-flexible/input/testTimePerTraversal/

Command:
` ./testTimePerTraversal.py [traversal=chosenTraversal] [tests=iterations of tests] [paths/to/yaml/files or/to/directories] `
Where traversal specifies a specific traversal and will run the scenario only for this traversal. This parameter may be left empty.

The tests parameter defines how often the scenario will be run per traversal. If it is not defined it will run once.

The last parameter needs to be the path to one or many input files.

### plotPerTraversalAndParameter.py
Plots all traversals from all input files according to the time and a defined parameter for the x-axis.
Traversals are colored according to their container.

Command: `./plotTuning.py parameter=[number, size, density, homogeneity] [path/To//Output/*.out ...]`

Wehere you can define the x-axis by setting the parameter. The parameter "number" will plot the time per number of particles,
"size" will plot per Boxsize, "density" per particle density and "homogeneity" the overall homogeneity of the scenario.
The last input parameter is the folder with one or many output files, produced by testTimePerTraversal.py

### outToCsvForPandas.py
Reads *.out files into csv which can by used for data analysis. Hardware has to be specified manually.

Command: `./outToCsvForPandas.py [hardware={name of hardware}] [path/To//Output/*.out ...]`


## Performance Measurements

### measurePerf.sh

### plotScript.gp

### strongscaling.sh
