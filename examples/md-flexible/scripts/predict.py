# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 31-07-2024

import pickle
import pandas as pd
import json
import os


def load_models_and_encoders(model_file: str) -> tuple:
    """
    Load the trained models and label encoders from a single pickle file.

    This function loads the previously saved RandomForest models and LabelEncoders from a single pickle file.
    These are used to make predictions based on live data.

    Args:
        model_file (str): The path to the pickle file containing the models and label encoders.

    Returns:
        models (dict): A dictionary containing the trained models for each target variable.
        label_encoders (dict): A dictionary containing LabelEncoders for decoding the predictions.
    """
    # Throw an error if the model file does not exist
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")

    # Load the models and label encoders from the pickle file
    with open(model_file, 'rb') as f:
        combined_data = pickle.load(f)
    # Return the models and label encoders
    return combined_data['models'], combined_data['label_encoders']


def preprocess_live_info(live_info: dict) -> pd.DataFrame:
    """
    Preprocess the live info received from C++ into a format suitable for model input.

    This function selects the relevant features from the live info dictionary (received from C++),
    formats them into a pandas DataFrame, and prepares them for input into the trained models.

    Args:
        live_info (dict): A dictionary containing the live info features from C++.

    Returns:
        pd.DataFrame: A DataFrame containing the processed features ready for prediction.
    """
    # Define the required features for the model input
    features = ['meanParticlesPerCell', 'maxParticlesPerCell', 'relativeParticlesPerCellStdDev', 'threadCount', 'relativeParticlesPerBlurredCellStdDev', 'skin']

    # Extract the relevant features from live info
    live_info_input = {feature: live_info[feature] for feature in features}

    # Convert the input into a pandas DataFrame with one row
    live_info_df = pd.DataFrame([live_info_input])

    return live_info_df


def predict(live_info: dict, models: dict, label_encoders: dict) -> dict:
    """
    Perform a forward pass through the trained models to predict the tuning configuration.

    This function takes the live info data, preprocesses it, and makes predictions using the loaded models.
    The predictions are then decoded back into their original labels using the corresponding LabelEncoders.

    Args:
        live_info (dict): A dictionary of live info data from C++.
        models (dict): A dictionary containing the trained models for each target variable.
        label_encoders (dict): A dictionary containing LabelEncoders for decoding predictions.

    Returns:
        dict: A dictionary containing the predicted tuning configuration (Container, Traversal, etc.).
    """
    live_info_df = preprocess_live_info(live_info)

    predictions_combined = {}

    # Perform prediction for each target
    for target, model in models.items():
        prediction_encoded = model.predict(live_info_df)
        prediction = label_encoders[target].inverse_transform(prediction_encoded)
        predictions_combined[target] = prediction[0]

    print(predictions_combined)

    predictions_split = predictions_combined['combined'].split(';')

    print(predictions_split)

    predictions = {
        'Container' : predictions_split[0],
        'Traversal' : predictions_split[1],
        'Load Estimator' : predictions_split[2],
        'Data Layout' : predictions_split[3],
        'Newton 3' : predictions_split[4],
        'CellSizeFactor' : predictions_split[5]
    }

    return predictions


def main(model_file: str, live_info_json: str) -> str:
    """
    Main execution function for performing a forward pass using live info.

    This function takes the model file name and live info JSON string as input, loads the models and encoders from the
    model file, and makes predictions for the tuning configuration based on the live info.

    Args:
        model_file (str): The file path of the saved models and label encoders.
        live_info_json (str): A JSON string representing the live info data.
    """
    # Load models and encoders
    models, label_encoders = load_models_and_encoders(model_file)

    # Parse the live info
    live_info = json.loads(live_info_json)

    # Make predictions
    predictions = predict(live_info, models, label_encoders)

    return json.dumps(predictions)
