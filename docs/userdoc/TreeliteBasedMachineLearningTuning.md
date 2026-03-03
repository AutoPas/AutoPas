# Treelite-based Machine Learning Tuning

AutoPas includes a Treelite-based machine learning tuning strategy: `treelite-based-decision-tree-tuning`.
This document describes how to build AutoPas for this strategy, generate training data, train and export the model, and use the generated Treelite artifacts at runtime.
The guide is tailored to HPC systems.

The central idea of the Treelite-based strategy is:

- train the decision-tree model in Python using `scikit-learn`
- export the trained model to Treelite
- run inference directly in C++ through Treelite's General Tree Inference Library (GTIL)

Compared to the Python-based tuning strategy, this avoids calling Python during prediction.

## Set-up

### Python

Python is still required for training and exporting the model, even though the runtime predictor itself is implemented in C++.

On HPC systems, this may mean loading a suitable Python module first.

### Python Virtual Environment and Dependencies

We recommend using a Python virtual environment to handle dependencies.
You can find details here: <https://docs.python.org/3/library/venv.html>

To create the environment:

```bash
python -m venv </path/to/new/virtual/environment>
```

To activate the environment:

```bash
source </path/to/new/virtual/environment>/bin/activate
```

Install the required Python packages:

```bash
pip install -r <autopas_src_dir>/src/autopas/tuning/tuningStrategy/decisionTreeTuning/requirements.txt
```

These dependencies are needed for:

- reading and merging AutoPas CSV log files
- training the `scikit-learn` model
- exporting the trained model to Treelite

## Generate the data

### AutoPas build

For the training workflow, you must enable the Live Info and Tuning Results loggers:

```bash
cmake .. -DAUTOPAS_LOG_LIVEINFO=ON -DAUTOPAS_LOG_TUNINGRESULTS=ON
```

Be careful when leaving the Live Info logger enabled outside of data generation, as it can noticeably reduce performance.

### Training data scenarios

How to generate the data is application-dependent.
It can be difficult to obtain a broad enough range of realistic training samples from real simulations alone, as a large range of data is required for successful training. Therefore, [Newcome et al., 2025](https://doi.org/10.1007/978-3-031-97635-3_35) proposed generating artificial particle distributions that are physically unrealistic but still represent a wide range of computationally realistic scenarios.

We provide an example for generating data for `md-flexible` in:
[`examples/md-flexible/input/trainingDataGeneration`](https://github.com/AutoPas/AutoPas/tree/master/examples/md-flexible/input/trainingDataGeneration)

See that directory's `README.md` for details.

## Training the model

The decision tree model is trained with the Python script:

[`src/autopas/tuning/tuningStrategy/decisionTreeTuning/train.py`](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/tuning/tuningStrategy/decisionTreeTuning/train.py)

Example:

```bash
python train.py --results-dir=/path/to/results/ --model-prefix=my_model
```

`/path/to/results/` should point to a root directory that contains the generated CSV files.
The script walks all leaf directories, matches Live Info and Tuning Results CSV files by rank and timestamp, and trains separate models for pairwise and triwise interactions when matching data is available.

By default, the script writes both Python and Treelite outputs, but only the Treelite outputs are needed for the Treelite-based tuning strategy.
The output targets can be controlled with:

- `--save-python` / `--no-save-python`
- `--save-treelite` / `--no-save-treelite`

### Output files

If Treelite export is enabled, the script may generate:

- `<model-prefix>_pairwise.tl`
- `<model-prefix>_pairwise_classes.txt`
- `<model-prefix>_triwise.tl`
- `<model-prefix>_triwise_classes.txt`
- `<model-prefix>_features.json`

The files are used as follows:

- `.tl`: serialized Treelite model checkpoint
- `_classes.txt`: one configuration class string per model output class
- `_features.json`: feature names in the exact order expected by the model

The Treelite export for an interaction type is skipped if:

- no training data was found for that interaction type
- the trained model has fewer than two classes

In that case, neither the Treelite model file nor the corresponding classes file is written for that interaction type.
The shared features file is only written if at least one Treelite model is exported.

## Using the tuning strategy

The following files are required at runtime for the Treelite-based tuning strategy:

- the Treelite model checkpoint (`.tl`)
- the corresponding classes file (`_classes.txt`)
- the shared features file (`_features.json`)

Set the `treelite-based-decision-tree-tuning` strategy in your AutoPas configuration.
This strategy is not intended to be combined with other tuning strategies.

On HPC systems, Python is only needed for the offline training and export step.
Once the Treelite artifacts have been generated, runtime prediction no longer depends on Python.
