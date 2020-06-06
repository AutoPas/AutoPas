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

## Compare Traversals per time

Enables you to compare several traversals per time and number of particles or domainsize or density.
In the plot you will find the time as the x-axis and the other parameter as y-axis.

Requirements:
* python3 (tested with 3.8.2)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.7.1)

### testTimePerTraversal.py
Examples can be found in AutoPas/examples/md-flexible/input/testTimePerTraversal/

Command:
` ./testTimePerTraversal.py [traversal=chosenTraversal] [tests=iterations of tests] [paths/to/yaml/files or/to/directories] `
Where traversal specifies a specific traversal and will run the scenario only for this traversal. This parameter may be left empty.

The tests parameter defines how often the scenario will be run per traversal. If it is not defined it will run once.

The last parameter needs to be the path to one or many input files.

### plotPerTraversalAndParameter.py
Plots all traversals from all input files according to the time and a defined parameter for the x-axis.
Traversals are colored according to their container.

Command: `./plotTuning.py parameter=[number, size, density] [path/To//Output/*.out ...]`

Wehere you can define the x-axis by setting the parameter. The parameter "number" will plot the time per number of particles,
"size" will plot per Boxsize and "density" per particle density.
The last input parameter is the folder with one or many output files, produced by testTimePerTraversal.py

## Performance Measurements

### measurePerf.sh

### plotScript.gp

### strongscaling.sh
