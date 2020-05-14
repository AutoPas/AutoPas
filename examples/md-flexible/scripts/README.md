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

Creates a plot of the predictions for each configuration for each md-flex output file given.

### plotDiffPredictionTest.py

Requirements:
* python3 (tested with Python 3.8.3rc1)
* [plotly](https://github.com/plotly/plotly.py) (tested with 4.7.1)

Creates a plot of the difference between the predictions and the tested times for each configuration for each md-flex output file given.

## Performance Measurements

### measurePerf.sh

### plotScript.gp

### strongscaling.sh
