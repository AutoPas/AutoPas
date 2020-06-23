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

## Performance Measurements

### measurePerf.sh

### plotScript.gp

### strongscaling.sh
