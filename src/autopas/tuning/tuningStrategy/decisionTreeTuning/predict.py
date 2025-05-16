# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 31-07-2024

import pickle
import pandas as pd
import json
import os
import numpy as np


def load_models_and_encoders(model_file: str) -> tuple:
    """
    Load the trained models, label encoders, and features from a single pickle file.

    This function loads the previously saved RandomForest models, LabelEncoders, and feature list
    from a single pickle file. These are used to make predictions based on live data.

    Args:
        model_file (str): The path to the pickle file containing the models, label encoders, and features.

    Returns:
        model (MultiOutputClassifier): The trained model for making predictions.
        label_encoders (dict): A dictionary containing LabelEncoders for decoding the predictions.
        features (list): A list of feature column names used during training.
    """
    # Throw an error if the model file does not exist
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")

    # Load the model, encoders, and features from the pickle file
    with open(model_file, 'rb') as f:
        combined_data = pickle.load(f)

    return combined_data['model'], combined_data['label_encoders'], combined_data['features']


def preprocess_live_info(live_info: dict, features: list) -> pd.DataFrame:
    """
    Preprocess the live info received from C++ into a format suitable for model input.

    This function selects the relevant features from the live info dictionary (received from C++),
    formats them into a pandas DataFrame, and prepares them for input into the trained models.

    Args:
        live_info (dict): A dictionary containing the live info features from C++.
        features (list): A list of feature column names to select from the live info.

    Returns:
        pd.DataFrame: A DataFrame containing the processed features ready for prediction.
    """
    # Extract the relevant features from live info
    live_info_input = {feature: live_info.get(feature, 0.0) for feature in features}

    # Convert the input into a pandas DataFrame with one row
    live_info_df = pd.DataFrame([live_info_input])

    return live_info_df


def predict(live_info: dict, model, label_encoders: dict, features: list) -> dict:
    """
    Perform a forward pass through the trained models to predict the tuning configuration.

    This function takes the live info data, preprocesses it, and makes predictions using the loaded models.
    The predictions are then decoded back into their original labels using the corresponding LabelEncoders.

    Args:
        live_info (dict): A dictionary of live info data from C++.
        model (MultiOutputClassifier): The trained model for making predictions.
        label_encoders (dict): A dictionary containing LabelEncoders for decoding predictions.
        features (list): A list of feature column names used for preprocessing the live info.

    Returns:
        dict: A dictionary containing the predicted tuning configuration (Container, Traversal, etc.) and confidence.
    """
    # Preprocess the live info into a DataFrame
    live_info_df = preprocess_live_info(live_info, features)

    # Get encoded predictions and probabilities
    prediction_encoded = model.predict(live_info_df)
    prediction_proba = model.predict_proba(live_info_df)

    predictions = {}
    total_confidence = 0
    num_targets = len(label_encoders)

    for i, target in enumerate(label_encoders.keys()):
        prediction = label_encoders[target].inverse_transform([prediction_encoded[0, i]])

        # Get the confidence score for the predicted label
        confidence_score = np.max(prediction_proba[i][0])
        total_confidence += confidence_score

        predictions[target] = prediction[0]

    # Calculate average confidence
    predictions["confidence"] = round(total_confidence / num_targets, 2)

    return predictions


def main(model_file: str, live_info_json: str) -> str:
    """
    Main execution function for performing a forward pass using live info.

    This function takes the model file name and live info JSON string as input, loads the models,
    encoders, and features from the model file, and makes predictions for the tuning configuration
    based on the live info.

    Args:
        model_file (str): The file path of the saved models, label encoders, and features.
        live_info_json (str): A JSON string representing the live info data.

    Returns:
        str: A JSON string representing the predicted tuning configuration and confidence score.
    """
    # Load models, encoders, and features
    model, label_encoders, features = load_models_and_encoders(model_file)

    # Parse the live info
    live_info = json.loads(live_info_json)

    # Make predictions
    predictions = predict(live_info, model, label_encoders, features)

    return json.dumps(predictions)
