#!/usr/bin/env python3
"""
Generates dummy pickle model files for DecisionTreeTuning unit tests.

Run automatically by CMake (generate_dt_test_models target) when
AUTOPAS_ENABLE_PYTHON_BASED_TUNING=ON.

Usage: python3 generate_test_models.py <output_dir>

Outputs
-------
test_model.pkl
    Valid model.  predict() returns the configuration the C++ tests assert:
    LinkedCells / lc_c08 / none / SoA / enabled / 1.0 / 1xVectorLength, with confidence 1.0.

test_model_invalid.pkl
    Structurally valid model (passes predict.py __init__) but whose label
    encoder produces a string with no semicolons, so predict() raises a
    ValueError when the number of split components does not match
    config_columns.  This exercises the error path in
    DecisionTreeTuning::getPredictionFromPython.
"""

import argparse
import os
import pickle

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

# Feature order must match liveInfoMap in DecisionTreeTuningTest.cpp.
FEATURES = [
    "avgParticlesPerCell",
    "maxParticlesPerCell",
    "homogeneity",
    "maxDensity",
    "particlesPerCellStdDev",
    "threadCount",
]

# Label components in the same order that train.py uses to build the alg_config column
# (Container;Traversal;LoadEstimator;DataLayout;Newton3;CellSizeFactor;VectorizationPattern).  String values come from
# the AutoPas Option::toString() maps in the C++ source.  If a new Configuration component is added to
# training/prediction, add it here and to CONFIG_COLUMNS in train.py.
_VALID_LABEL_COMPONENTS = [
    ("Container",      "LinkedCells"),
    ("Traversal",      "lc_c08"),
    ("Load Estimator", "none"),
    ("Data Layout",    "SoA"),
    ("Newton 3",       "enabled"),
    ("CellSizeFactor", "1.0"),
    ("VectorizationPattern", "1xVectorLength")
]
VALID_LABEL = ";".join(v for _, v in _VALID_LABEL_COMPONENTS)

# The ordered configuration component names stored alongside the model, mirroring CONFIG_COLUMNS in train.py. predict.py
# reads this to split a predicted label back into its components, so it must be present in every model file.
CONFIG_COLUMNS = [name for name, _ in _VALID_LABEL_COMPONENTS]

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


def _dump(path: str, pairwise_model, pairwise_label_encoder) -> None:
    data = {
        "pairwise_model": pairwise_model,
        "pairwise_label_encoder": pairwise_label_encoder,
        "triwise_model": None,
        "triwise_label_encoder": None,
        "features": FEATURES,
        "config_columns": CONFIG_COLUMNS,
    }
    with open(path, "wb") as fh:
        pickle.dump(data, fh)


def create_valid_model(output_dir: str) -> None:
    """Predict VALID_LABEL with confidence 1.0 for any input."""
    model, le = _make_model(VALID_LABEL)
    _dump(os.path.join(output_dir, "test_model.pkl"), model, le)


def create_invalid_model(output_dir: str) -> None:
    """
    Passes predict.py __init__ but causes a ValueError in predict() because the
    decoded label has no semicolons, so the number of split components (1) does
    not match config_columns.
    """
    model, le = _make_model("invalid_label_no_semicolons")
    _dump(os.path.join(output_dir, "test_model_invalid.pkl"), model, le)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", help="Directory to write the .pkl files into")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    create_valid_model(args.output_dir)
    create_invalid_model(args.output_dir)
    print(f"Test models written to {args.output_dir}")


if __name__ == "__main__":
    main()