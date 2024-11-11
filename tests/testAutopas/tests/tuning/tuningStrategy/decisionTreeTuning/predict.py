# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 20-09-2024

import json
import sys


def main(model_file: str, live_info_json: str) -> str:
    """
    Main function that receives the model file and live_info as arguments, simulates a mock prediction,
    and returns the predicted configuration.

    Args:
        model_file (str): Path to the model file (not used in this mock).
        live_info_json (str): JSON string containing the live info from the C++ side.

    Returns:
        str: JSON string representing the predicted configuration.
    """
    if model_file == "test_model.pkl":
        predicted_config = {
            "Container": "LinkedCells",
            "Traversal": "lc_c08",
            "Load Estimator": "none",
            "Data Layout": "SoA",
            "Newton 3": "enabled",
            "CellSizeFactor": "1.0"
        }
    elif model_file == "test_model_invalid_type.pkl":
        predicted_config = {
            "Container": "LinkedCells",
            "Traversal": "lc_c08",
            "Load Estimator": "none",
            "Data Layout": "SoA",
            "Newton 3": "enabled",
            "CellSizeFactor": "invalid_number"
        }
    else:
        raise ValueError("Invalid model file")
    # Convert the predicted configuration to a JSON string and return it
    return json.dumps(predicted_config)
