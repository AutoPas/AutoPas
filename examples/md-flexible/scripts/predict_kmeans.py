import json

import joblib
import pandas as pd
import os

def load_model_and_cluster_configuration(model_path: str, cluster_configuration_path: str) -> tuple:
    """
    This function loads the k-means prediction pipeline from the model_path and the configuration-cluster mapping from the path in
    cluster_configuration_path

    :param model_path:
    :param cluster_configuration_path:
    :return:
    """
    #look if the path exists for the model and the configurations mapping
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    if not os.path.exists(cluster_configuration_path):
        raise FileNotFoundError(f"Cluster_configuration file not found: {cluster_configuration_path}")

    #load the pipeline cluster model
    model = joblib.load(model_path)

    #open the configuration mapping
    with open(cluster_configuration_path, 'r') as f:
        config_mapping = json.load(f)

    #json keys are strings cast them to int
    config_mapping = {int(k): v for k, v in config_mapping.items()}

    return model, config_mapping



def preprocess_live_info(live_info: dict) -> pd.DataFrame:
    """
    This function preprocesses the liveInfo by extracting the relevant features.
    :param live_info:
    :return:
    """

    #the relevant columns that are used to predict the configuration
    features = ['minParticlesPerCell', 'meanParticlesPerCell', 'maxParticlesPerCell', 'relativeParticlesPerCellStdDev',
                                 'relativeParticlesPerBlurredCellStdDev', 'numOwnedParticles', 'medianParticlesPerCell']

    #read in the values in the feature columns in case of an key_error return 0.0 as default value
    data = {feature: live_info.get(feature, 0.0) for feature in features}

    live_info_df = pd.DataFrame([data])

    #reindex the columns to make sure they are in correct order
    return live_info_df.reindex(columns=features)



def predict(live_info: dict, pipeline, config_mapping) -> str:
    """
    This function uses the k-means pipeline to predict a cluster and then returns a configuration that is used in that cluster.

    :param live_info:
    :param pipeline:
    :param config_mapping:
    :return:
    """
    #preprocess the liveInfo
    live_info_input = preprocess_live_info(live_info)

    #the model predicts the cluster based on the liveInfo
    cluster_label = pipeline.predict(live_info_input)[0]

    #the configurations that are in that cluster are chosen
    configurations = config_mapping[cluster_label]

    return configurations



def main(model_path: str, cluster_configuration_path: str, live_info_json: str) -> str:
    """
    This function is called by the cluster-based-tuning strategy it load the prediction model and the configuration mapping
    and then makes a prediction based on the live_info_json provided as an argument.

    :param model_path:
    :param cluster_configuration_path:
    :param live_info_json:
    :return:
    """
    #loads the model and the configuration mapping
    model, config_mapping = load_model_and_cluster_configuration(model_path, cluster_configuration_path)

    #convert the liveInfo into json format
    live_info = json.loads(live_info_json)

    #makes a prediction and returns suitable configurations
    configuration_predictions = predict(live_info, model, config_mapping)

    #serialize to json
    configuration_predictions_json = json.dumps(configuration_predictions)

    return configuration_predictions_json
