# Python-based Machine Learning Tuning

AutoPas includes a python-based machine learning tuning strategy: decision-tree-tuning. This document provides a guide on how to set this up, generate the data, training the model, and use the tuning strategy. The guide is tailored to HPC systems.

The central premise of python-based tuning strategies is to use PyBind11 is call python from within AutoPas. That allows for easier experimentation and use of Python machine learning libraries such as sklearn.

We, therefore, will also discuss possibilities to modify and experiment with the python-based machine learning tuning strategy beyond the random-forest based strategy.

## Set-up

### Python

To use this tuning strategy, you will need a suitable version of Python, compatible with the version of PyBind11 being used. This means, currently, at least v3.7. (See [PyBind11 releases](https://github.com/pybind/pybind11/releases) and cross-reference against the version bundled under `./libs` or your own version for the most-up-to-date information.) On HPC systems, this sometimes means the loading a module for a more up-to-date version of Python.

We have tested this with Python 3.11.7.

### Python Virtual Enviroment and Dependencies

We recommend using a Python virtual enviroment to handle dependencies. You can find details of this [here](https://docs.python.org/3/library/venv.html).

To create the enviroment:
```bash
python -m venv </path/to/new/virtual/environment>
```

To activate the enviroment:
```bash
source </path/to/new/virtual/environment>/bin/activate
```

Install dependencies using pip. For tested dependencies required by the decision tree tuning strategy, you can install these automatically:

```bash
pip install -r <autopas_src_dir>/src/autopas/tuning/tuningStrategy/decisionTreeTuning/requirements.txt
```

## Generate the data

### AutoPas build

For the decision tree tuning strategy, you must enable the live info and tuning results loggers (in the build directory):

```bash
cmake .. -DAUTOPAS_LOG_LIVEINFO=ON -DAUTOPAS_LOG_TUNINGRESULTS=ON
```

Be careful of keeping the live info logger enabled beyond data generation, as it can have a significant negative performance impact.

### Training data scenarios

How to generate the data is ultimately application dependent. It is potentially challenging to generate appropriate training data through real simulation, as a large range of data is required for successful training. Therefore, [Newcome et al., 2025](https://doi.org/10.1007/978-3-031-97635-3_35) proposed creating fake simulation data by artifically distributing particles in such a way that, whilst physically inaccurate, recreated a wide range of computationally realistic scenarios (e.g. realistic numbers of particles per cell, but some unrealistically close together - which makes no impact on algorithmic performance). 

We provide an example for generating data for the md-flexible example code. See the directory [`examples/md-flexible/input/trainingDataGeneration`](https://github.com/AutoPas/AutoPas/tree/master/examples/md-flexible/input/trainingDataGeneration) and in particular its `README.md`, where more details are given.

Newcome et al., 2025 pre-dated [dynamically-rebuilt containers](https://github.com/AutoPas/AutoPas/pull/821), and so the effects of variable rebuild frequencies was not considered. Therefore, the success of the above methodology may be limited, and will be addressed in the future.

## Training the model

The decision tree tuning model can be trained with the python script [`src/autopas/tuning/tuningStrategy/decisionTreeTuning/train.py`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/tuning/tuningStrategy/decisionTreeTuning/train.py).

For example

```bash
python train.py --results-dir=/path/to/results/
```

`/path/to/results/` should be a single directory containing all result files. Live Info and Tuning Result logs will be matched together based on their shared time-and-rank-stamps. As it is possible for multiple experiments to run with the same time stamp, we therefore recommend running a single experiment per root directory, such as with the example training data generators.

## Using the tuning strategy

Set the decision-tree-tuning strategy. This is not designed to be used with any other tuning strategies.

If using a HPC system, remember to, in the jobscript file, load the up-to-date version of python and activate the virtual enviroment.



