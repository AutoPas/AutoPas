import json

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import joblib
import os
import glob
from sklearn.metrics import silhouette_score, make_scorer


def load_data_from_directory_drop(results_dir: str) -> pd.DataFrame:
    """
    Load simulation data from the results_dir into one DataFrame-

    This function recursively iterates trough the directory that is passed as an argument.
    The files stored there should be created by the LiveInfoLogger and the TuningResultsLogger and be of type csv.
    They are then merged and duplicates with the same iteration number are removed.
    :param results_dir:
    :return:
    """
    merged_dfs = []

    # Loop through all leaf folders i.e. the ones in which the csv files are outputted in

    for root, dirs, files in os.walk(results_dir):
        if not dirs:
            # Get Live Info file

            live_info_file_pattern = os.path.join(root, "*liveInfoLogger*")
            live_info_file_list = glob.glob(live_info_file_pattern)

            if len(live_info_file_list) == 0:
                print(f"No Live Info csv file found in folder {root}")

            if len(live_info_file_list) > 1:
                print(f"More than one Live Info csv file found in folder {root}")

            live_info_file = live_info_file_list[0]

            # Get Iteration file

            tuning_results_file_pattern = os.path.join(root, "*tuningResults*")
            tuning_results_file_list = glob.glob(tuning_results_file_pattern)

            if len(tuning_results_file_list) == 0:
                print(f"No Iteration csv file found in folder {root}")

            if len(tuning_results_file_list) > 1:
                print(f"More than one Iteration csv file found in folder {root}")

            tuning_results_file = tuning_results_file_list[0]

            # Load the live info and tuning results
            live_info_df = pd.read_csv(live_info_file)

            # remove duplicates with the same iteration number
            live_info_df = live_info_df.drop_duplicates(subset='Iteration', keep='first')

            tuning_results_df = pd.read_csv(tuning_results_file)

            # remove duplicates with the same iteration number
            tuning_results_df = tuning_results_df.drop_duplicates(subset='Iteration', keep='first')

            # Merge them on 'Iteration' column
            merged_df = pd.merge(live_info_df, tuning_results_df, on='Iteration', how='right')

            # Clean column names
            merged_df.columns = merged_df.columns.str.strip()

            # Append the merged DataFrame to the list
            merged_dfs.append(merged_df)

    # Combine all merged DataFrames into one
    final_merged_df = pd.concat(merged_dfs, ignore_index=True)

    return final_merged_df


def preprocess_data(merged_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    This function is used to preprocess the training data, it selects the training features and creates a training dataset,
    and a combined dataset that contains the liveInfo columns for training and the
    configurations that were used with that LiveInfo in the simulation.

    :param merged_df:
    :return:
    """
    # the features used for training k-means (X_train)
    train_features = ['minParticlesPerCell', 'meanParticlesPerCell', 'maxParticlesPerCell',
                      'relativeParticlesPerCellStdDev',
                      'relativeParticlesPerBlurredCellStdDev', 'numOwnedParticles', 'medianParticlesPerCell']

    # tuning features
    tuning_features = ['Traversal', 'Newton 3', 'Container', 'Data Layout']

    # this is used for a dataframe that contains both the liveInfo and the corresponding tuning results
    combined_features = train_features + tuning_features

    # X_train used to train k-means
    X_train = merged_df[train_features].copy()

    combined_df = merged_df[combined_features].copy()

    return X_train, combined_df


def train_model(X_train: pd.DataFrame, combined_df: pd.DataFrame) -> tuple:
    """
    Train a k-means pipeline with a StandardScaler preprocessor and remove outliers with PCA.

    This function checks how many different configurations of there are, and then trains a k-means pipeline with that number of unique clusters.

    :param X_train:
    :param combined_df:
    :return:
    """
    # check how many different tuning configurations there are
    num_configurations = combined_df[['Traversal', 'Container', 'Newton 3', 'Data Layout']].apply(tuple,
                                                                                                  axis=1).nunique()

    # create a pipeline to standardize the data, remove outliers with pca and group with k-means
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA(n_components=0.95)),  # optional
        ('kmeans', KMeans(n_clusters=num_configurations,
                          init='k-means++',
                          n_init=10,
                          max_iter=300,
                          random_state=42))
    ])

    # fit the pipeline to the liveInfo data that is contained in X_train
    pipeline.fit(X_train)

    # k-means predicts and groups the X_train data into groups that should have similar configurations
    combined_df['configuration_group'] = pipeline.predict(X_train)

    return pipeline, combined_df


def save_model_and_configurations(pipeline, combined_df: pd.DataFrame):
    """
    This function saves the pipeline and creates a dict, for each cluster a list of all unique configurations that are assigned to that cluster is saved.
    This creates the configuration mapping.

    :param pipeline:
    :param combined_df:
    :return:
    """
    # the relevant columns
    config_cols = ['Traversal', 'Container', 'Newton 3', 'Data Layout']

    # all unique cluster groups
    groups = combined_df['configuration_group'].unique()

    # dict: cluster_id -> list of all unique configurations in that cluster
    config_mapping = {
        int(g): combined_df[combined_df['configuration_group'] == g]
        [config_cols]
        .drop_duplicates()
        .to_dict(orient='records')
        for g in groups
    }

    joblib.dump(pipeline, 'cluster-model.pkl')

    with open('cluster_configuration_mapping.json', 'w') as f:
        json.dump(config_mapping, f, indent=2)

if __name__ == '__main__':
    """
    The main function loads LiveInfo and Tuning Results Data from the data directory in the project root, then it 
    preprocesses and trains a pipeline with k-means and saves the model and the corresponding configurations
    to the examples/md-flexible/scripts directory
    """

    #load the simulation data
    merged_df = load_data_from_directory_drop('../../../data')

    #preprocess data
    X_train, combined_df = preprocess_data(merged_df)

    #train the model
    pipeline, combined_df = train_model(X_train, combined_df)

    #save model and configuration mapping
    save_model_and_configurations(pipeline, combined_df)
