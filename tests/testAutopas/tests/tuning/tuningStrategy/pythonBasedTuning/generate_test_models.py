#!/usr/bin/env python3
"""
Generates dummy pickle model files for PythonBasedTuning unit tests.

Outputs
-------
test_model.pkl
    Valid model.  predict() returns only the configuration whose values the C++ test passed.

test_model_invalid.pkl
    Structurally valid model (passes predict.py __init__) but whose label
    encoder produces a string with no semicolons, so predict() raises a
    ValueError when the number of split components does not match
    config_columns.  This exercises the error path in
    PythonBasedTuning::getPredictionFromPython.

test_model_ranked.pkl
    Valid model.  predict() returns the configurations whose values the C++
    test passed, in the given order, and leaves out one with probability zero.
"""

import os
import pickle

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

# Feature order must match liveInfoMap in PythonBasedTuningTest.cpp.
FEATURES = [
    "avgParticlesPerCell",
    "maxParticlesPerCell",
    "homogeneity",
    "maxDensity",
    "particlesPerCellStdDev",
    "threadCount",
]

# Single training sample – the same values the test supplies as live info.
TRAINING_INPUT = np.array([[6.82, 33.0, 0.42, 1.17, 0.03, 1.0]])


def _make_model(label: str):
    """Return a (RandomForestClassifier, LabelEncoder) pair trained on one sample."""
    le = LabelEncoder()
    le.fit([label])
    y = le.transform([label])
    clf = RandomForestClassifier(n_estimators=10, random_state=42)
    clf.fit(TRAINING_INPUT, y)
    return clf, le


def _dump(path: str, pairwise_model, pairwise_label_encoder, config_columns: list) -> None:
    data = {
        "pairwise_model": pairwise_model,
        "pairwise_label_encoder": pairwise_label_encoder,
        "triwise_model": None,
        "triwise_label_encoder": None,
        "features": FEATURES,
        "config_columns": config_columns,
    }
    with open(path, "wb") as fh:
        pickle.dump(data, fh)


def create_valid_model(output_dir: str, config_columns: list, config_values: list) -> None:
    """
    Predict only the configuration with config_values for any input. The label joins the values with ';', like
    build_alg_config in train.py.
    """
    model, le = _make_model(";".join(config_values))
    _dump(os.path.join(output_dir, "test_model.pkl"), model, le, config_columns)


def create_ranked_model(output_dir: str, config_columns: list, ranked_config_values: list,
                        zero_probability_config_values: list) -> None:
    """
    For the live info TRAINING_INPUT, predict the configurations with ranked_config_values from most to least probable,
    and the configuration with zero_probability_config_values with probability zero.

    Without bootstrapping, every tree is fit on all samples, and the identical samples at TRAINING_INPUT end up in one
    leaf. So the probabilities there are exactly the label frequencies among those samples.
    """
    ranked_labels = [";".join(values) for values in ranked_config_values]
    # Each label gets one sample less than the label before it, so the probabilities strictly decrease.
    labels = [label for rank, label in enumerate(ranked_labels) for _ in range(len(ranked_labels) - rank)]
    inputs = np.repeat(TRAINING_INPUT, len(labels), axis=0)
    # The zero-probability configuration is only trained far from TRAINING_INPUT.
    labels.append(";".join(zero_probability_config_values))
    inputs = np.vstack([inputs, TRAINING_INPUT * 10])

    label_encoder = LabelEncoder()
    encoded_labels = label_encoder.fit_transform(labels)
    model = RandomForestClassifier(n_estimators=10, random_state=42, bootstrap=False)
    model.fit(inputs, encoded_labels)
    _dump(os.path.join(output_dir, "test_model_ranked.pkl"), model, label_encoder, config_columns)


def create_invalid_model(output_dir: str, config_columns: list) -> None:
    """
    Passes predict.py __init__ but causes a ValueError in predict() because the
    decoded label has no semicolons, so the number of split components (1) does
    not match config_columns.
    """
    model, le = _make_model("invalid_label_no_semicolons")
    _dump(os.path.join(output_dir, "test_model_invalid.pkl"), model, le, config_columns)