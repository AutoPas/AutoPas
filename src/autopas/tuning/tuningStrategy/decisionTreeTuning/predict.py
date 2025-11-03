# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 31-07-2024

import pickle
import pandas as pd
import json
import os
import numpy as np


def load_model_and_encoder(model_file: str, interaction_type: str) -> tuple:
    """
    Load the trained model, label encoder, and features from a single pickle file.

    This function loads the previously saved RandomForest model, LabelEncoder, and feature list
    from a single pickle file. These are used to make predictions based on live data.

    Args:
        model_file (str): The path to the pickle file containing the model, label encoder, and features.
        interaction_type (str): The interaction type of the algorithm to be predicted.

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
        try:
            combined_data = pickle.load(f)
        except Exception as e:
            raise Exception(f"Error loading pickle file {model_file}: {e}")

    if interaction_type == "pairwise":
        try:
            model = combined_data['pairwise_model']
            label_encoder = combined_data['pairwise_label_encoder']
        except:
            raise Exception("The Pairwise Model or Label Encoder was not found.")


    elif interaction_type == "triwise":
        try:
            model = combined_data['triwise_model']
            label_encoder = combined_data['triwise_label_encoder']
        except:
            raise Exception("The Triwise Model or Label Encoder was not found.")

    else:
        raise ValueError(f"Invalid interaction type: {interaction_type}")

    try:
        features = combined_data['features']
    except:
        raise Exception("The Features were not found.")

    # Check that the model exists (e.g. that a triwise model was actually trained)
    if model is None:
        raise Exception(f"There is no model for {interaction_type} interactions. This could be because no Tuning Results"
                        f" logs for this interaction type were used during training.")

    return model, label_encoder, features


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


def predict(live_info: dict, model, label_encoder, features: list) -> dict:
    """
    Perform a forward pass through the trained models to predict the tuning configuration.

    This function takes the live info data, preprocesses it, and makes predictions using the loaded model.
    The prediction is then decoded back into its original label using the corresponding LabelEncoder, and returned
    together with the model's "confidence"/"probability". (The leaves of each decision tree contain a number of training
    samples. The "probability" of each configuration is simply the percentage of samples in that leaf where the
    configuration was chosen. The random forest simply averages these probabilities.)

    Args:
        live_info (dict): A dictionary of live info data from C++.
        model (RandomForestClassifier): The trained model for making predictions.
        label_encoder (LabelEncoder): A LabelEncoder for decoding predictions.
        features (list): A list of feature column names used for preprocessing the live info.

    Returns:
        prediction (dict): Dictionary containing the predicted optimal algorithmic configuration and the confidence
        score.
    """
    # Preprocess the live info into a DataFrame
    live_info_df = preprocess_live_info(live_info, features)

    # Get prediction probabilities
    probabilities = model.predict_proba(live_info_df)[0]

    # Todo, we may be able to improve the generalizability of this tuning strategy by taking e.g. the top 5 predictions
    predicted_class_id = np.argmax(probabilities)

    prediction = label_encoder.inverse_transform([predicted_class_id])[0]
    prediction_split = prediction.split(";")

    prediction = {
        "Container": prediction_split[0],
        "Traversal": prediction_split[1],
        "Load Estimator": prediction_split[2],
        "Data Layout": prediction_split[3],
        "Newton 3": prediction_split[4],
        "CellSizeFactor": prediction_split[5],
        "confidence": probabilities[predicted_class_id]
    }

    return prediction


def main(model_file: str, live_info_json: str, interaction_type: str) -> str:
    """
    Main execution function for performing a forward pass using live info.

    This function takes the model file name and live info JSON string as input, loads the model,
    encoder, and features from the model file, and makes predictions for the tuning configuration
    based on the live info.

    Args:
        model_file (str): The file path of the saved models, label encoders, and features.
        live_info_json (str): A JSON string representing the live info data.
        interaction_type (str): The interaction type of the algorithm to be predicted.

    Returns:
        str: A JSON string representing the predicted tuning configuration and confidence score.
    """
    # Load model, encoder, and features
    model, label_encoder, features = load_model_and_encoder(model_file, interaction_type)

    # Parse the live info
    live_info = json.loads(live_info_json)

    # Make prediction
    prediction = predict(live_info, model, label_encoder, features)

    return json.dumps(prediction)
