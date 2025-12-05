# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 31-07-2024

import pickle
import pandas as pd
import json
import os
import numpy as np

class DecisionTreeTuning:
    """
    Stores and predicts with a trained Decision Tree model, label encoder, and the live info features that it
    takes as input.

    This class is designed to mimic the C++ class, i.e. be initialized when that class is constructed, and the predict
    method used when getPredictionFromPython is called.
    """
    def __init__(self, model_file: str, interaction_type: str):
        """
        Load the trained model, label encoder, and features from a pickle file.

        Args:
            model_file (str): The path to the pickle file containing the model, label encoder, and features.
            interaction_type (str): The interaction type of the algorithm to be predicted.
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
                self.model = combined_data['pairwise_model']
                self.label_encoder = combined_data['pairwise_label_encoder']
            except:
                raise Exception("The Pairwise Model or Label Encoder was not found.")


        elif interaction_type == "triwise":
            try:
                self.model = combined_data['triwise_model']
                self.label_encoder = combined_data['triwise_label_encoder']
            except:
                raise Exception("The Triwise Model or Label Encoder was not found.")

        else:
            raise ValueError(f"Invalid interaction type: {interaction_type}")

        try:
            self.features = combined_data['features']
        except:
            raise Exception("The Features were not found.")

        # Check that the model exists (e.g. that a triwise model was actually trained)
        if self.model is None:
            raise Exception(f"There is no model for {interaction_type} interactions. This could be because no Tuning Results"
                            f" logs for this interaction type were used during training.")


    def preprocess_live_info(self, live_info_json: str) -> pd.DataFrame:
        """
        Preprocess the live info received from C++ into a format suitable for model input.

        This function selects the relevant features from the live info dictionary (received from C++),
        formats them into a pandas DataFrame, and prepares them for input into the trained models.

        Args:
            live_info_json (str): A JSON containing the live info features from C++.

        Returns:
            pd.DataFrame: A DataFrame containing the processed features ready for prediction.
        """
        # Convert the live info JSON string into a dictionary
        live_info = json.loads(live_info_json)

        # Extract the relevant features from live info
        live_info_input = {feature: live_info.get(feature, 0.0) for feature in self.features}

        # Convert the input into a pandas DataFrame with one row
        live_info_df = pd.DataFrame([live_info_input])

        return live_info_df


    def predict(self, live_info_json: str) -> str:
        """
        Perform a forward pass through the trained models to predict the tuning configuration.

        This function takes the live info data, preprocesses it, and makes predictions using the loaded model.
        The prediction is then decoded back into its original label using the corresponding LabelEncoder, and returned
        together with the model's "confidence"/"probability". (The leaves of each decision tree contain a number of training
        samples. The "probability" of each configuration is simply the percentage of samples in that leaf where the
        configuration was chosen. The random forest simply averages these probabilities.)

        Args:
            live_info_json (str): A JSON of live info data from C++.

        Returns:
            str: JSON containing the predicted optimal algorithmic configuration and confidence score.
        """
        # Preprocess the live info into a DataFrame
        live_info_df = self.preprocess_live_info(live_info_json)

        # Get prediction probabilities
        probabilities = self.model.predict_proba(live_info_df)[0]

        # Todo, we may be able to improve the generalizability of this tuning strategy by taking e.g. the top 5 predictions
        predicted_class_id = np.argmax(probabilities)

        prediction = self.label_encoder.inverse_transform([predicted_class_id])[0]
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

        return json.dumps(prediction)


