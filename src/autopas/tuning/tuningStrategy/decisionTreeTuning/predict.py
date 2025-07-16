# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 31-07-2024

import pickle
import pandas as pd
import json
import os
import numpy as np


def load_model_and_encoder(model_file: str) -> tuple:
    """
    Load the trained model, label encoder, and features from a single pickle file.

    This function loads the previously saved RandomForest model, LabelEncoder, and feature list
    from a single pickle file. These are used to make predictions based on live data.

    Args:
        model_file (str): The path to the pickle file containing the model, label encoder, and features.

    Returns:
        model (RandomForestClassifier): The trained model for making predictions.
        label_encoder (LabelEncoder): A LabelEncoder for decoding the predictions.
        features (list): A list of feature column names used during training.
    """
    # Throw an error if the model file does not exist
    if not os.path.exists(model_file):
        raise FileNotFoundError(f"Model file not found: {model_file}")

    # Load the model, encoders, and features from the pickle file
    with open(model_file, 'rb') as f:
        combined_data = pickle.load(f)

    return combined_data['model'], combined_data['label_encoder'], combined_data['features']


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


def predict(live_info: dict, model, label_encoder: dict, features: list) -> dict:
    """
    Perform a forward pass through the trained models to predict the tuning configuration.

    This function takes the live info data, preprocesses it, and makes predictions using the loaded models.
    The prediction is then decoded back into its original label using the corresponding LabelEncoder.

    Args:
        live_info (dict): A dictionary of live info data from C++.
        model (RandomForestClassifier): The trained model for making predictions.
        label_encoder (LabelEncoder): A LabelEncoder for decoding predictions.
        features (list): A list of feature column names used for preprocessing the live info.

    Returns:
        prediction (str): The predicted optimal algorithmic configuration.
        prediction_confidence (float): The probability/confidence that the input belongs to the predicted alg. conf.
        class.
    """
    # Preprocess the live info into a DataFrame
    live_info_df = preprocess_live_info(live_info, features)

    # Get encoded predictions and probabilities
    prediction_encoded = model.predict(live_info_df)

    prediction_combined = label_encoder.inverse_transform([prediction_encoded[0]])[0]
    prediction_split = prediction_combined.split(';')
    prediction = {
        'Container' : prediction_split[0],
        'Traversal' : prediction_split[1],
        'Load Estimator' : prediction_split[2],
        'Data Layout' : prediction_split[3],
        'Newton 3' : prediction_split[4],
        'CellSizeFactor' : prediction_split[5]
    }
    
    prediction["confidence"]  = np.max(model.predict_proba(live_info_df)[0])

    return prediction


def main(model_file: str, live_info_json: str) -> str:
    """
    Main execution function for performing a forward pass using live info.

    This function takes the model file name and live info JSON string as input, loads the model,
    encoder, and features from the model file, and makes predictions for the tuning configuration
    based on the live info.

    Args:
        model_file (str): The file path of the saved models, label encoders, and features.
        live_info_json (str): A JSON string representing the live info data.

    Returns:
        str: A JSON string representing the predicted tuning configuration and confidence score.
    """
    # Load model, encoder, and features
    model, label_encoder, features = load_model_and_encoder(model_file)

    # Parse the live info
    live_info = json.loads(live_info_json)

    # Make prediction
    prediction = predict(live_info, model, label_encoder, features)

    return json.dumps(prediction)
