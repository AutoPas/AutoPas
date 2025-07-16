# @File: train.py
# @Author Abdulkadir Pazar
# @Date 10-08-2024

import os
import argparse
import glob
from pyexpat import features

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from sklearn.multioutput import MultiOutputClassifier
import pickle
import numpy as np

def extract_identifier(filename: str) -> str:
    """
    Extract the rank and timestamp from the filename, to serve as a unique identifier to match Live Info csvs to
    tuning results csvs.

    I.e. 'AutoPas_liveInfoLogger_Rank0_2025-05-13_14-11-48.csv' becomes 'Rank0_2025-05-13_14-11-48.csv'

    Args:
        filename (str): The filename from which to extract the identifier.

    Returns:
        str: The extracted identifier.
    """
    return filename.split('_', 2)[2]


def load_data_from_directory(results_dir: str) -> pd.DataFrame:
    """
    Load and merge multiple live info and tuning results CSV files assuming leaf folders contain all experiments with
    each folder containing matching live info and tuning result CSV files. Files are matched together using their rank
    and timestamps.

    Args:
        results_dir (str): The directory containing all experiment CSV files.

    Returns:
        pd.DataFrame: A merged DataFrame containing the combined live info and tuning results from all files.
    """
    merged_dfs = []

    # Loop through all leaf folders i.e. the ones in which the csv files are outputted in

    for root, dirs, files in os.walk(results_dir):
        if not dirs:
            # Get Live Info files and sort them
            live_info_file_pattern = os.path.join(root,"*liveInfoLogger*")
            live_info_file_list = glob.glob(live_info_file_pattern)

            # Get Tuning Results files
            tuning_results_file_pattern = os.path.join(root,"*tuningResults*")
            tuning_results_file_list = glob.glob(tuning_results_file_pattern)

            # Create dictionaries mapping identifiers to files and use it to get a list of matching file pairs
            live_info_map = {extract_identifier(f): f for f in live_info_file_list}
            tuning_results_map = {extract_identifier(f): f for f in tuning_results_file_list}

            common_identifiers = set(live_info_map) & set(tuning_results_map)
            csv_file_pairs = [(live_info_map[id], tuning_results_map[id]) for id in common_identifiers]

            for [live_info_file, tuning_results_file] in csv_file_pairs:
                # Load the live info and tuning results
                live_info_df = pd.read_csv(live_info_file)
                tuning_results_df = pd.read_csv(tuning_results_file)

                # Combine the algorithmic configuration into one column
                tuning_results_df['alg_config'] = (tuning_results_df['Container'] + ';'
                                                   + tuning_results_df['Traversal']
                                                   + ';' + tuning_results_df['Load Estimator'] + ';'
                                                   + tuning_results_df['Data Layout'] + ';'
                                                   + tuning_results_df['Newton 3'] + ';'
                                                   + tuning_results_df['CellSizeFactor'].map(str))

                # Merge them on 'Iteration' column
                merged_df = pd.merge(live_info_df, tuning_results_df, on='Iteration', how='right')

                # Clean column names
                merged_df.columns = merged_df.columns.str.strip()

                # Append the merged DataFrame to the list
                merged_dfs.append(merged_df)

    # Combine all merged DataFrames into one
    final_merged_df = pd.concat(merged_dfs, ignore_index=True)

    return final_merged_df


def preprocess_data(merged_df: pd.DataFrame, features: list) -> tuple:
    """
    Preprocess the merged DataFrame for model training.

    This function selects the relevant features from the live info and the algorithm configuration.
    It encodes a categorical target variable using LabelEncoder and splits the merged DataFrame into features (X)
    and the target label (y).

    Args:
        merged_df (pd.DataFrame): The merged DataFrame containing live info and tuning results.
        features (list): List of column names to be used as features.

    Returns:
        X (pd.DataFrame): A DataFrame of feature columns used for training.
        y (pd.DataFrame): A DataFrame of the alg_config target column to be predicted.
        label_encoder (LabelEncoder): A LabelEncoder used for encoding the target variables.
    """
    label_encoder = LabelEncoder()

    merged_df['alg_config'] = label_encoder['alg_config'].fit_transform(merged_df['alg_config'])
    X = merged_df[features]
    y = merged_df['alg_config']
    return X, y, label_encoder


def train_model(X: pd.DataFrame, y: pd.DataFrame, test_size: float, n_estimators: int) -> RandomForestClassifier:
    """
    Train a RandomForestClassifier.

    This function splits the data into training and testing sets and trains a RandomForestClassifier on the training
    data. It evaluates the model by calculating the accuracy on the test data, but be cautious using this result, as
    an incorrect prediction does not matter much if its performance is not much worse than the optimum.

    Args:
        X (pd.DataFrame): The feature set for training the model.
        y (pd.DataFrame): The target set for training the model.
        test_size (float): Fraction of data to be used for testing.
        n_estimators (int): Number of estimators for the RandomForestClassifier.

    Returns:
        RandomForestClassifier: The trained RandomForestClassifier model.
    """
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=42)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    accuracy = accuracy_score(y_test.iloc[:, 'alg_config'], y_pred[:, 'alg_config'])
    print(f'Accuracy: {accuracy:.3f}')

    return model


def save_model_and_encoder(model: RandomForestClassifier, label_encoder: dict, features: list, output_file: str) -> None:
    """
    Save the trained model, LabelEncoder, and feature list to a file.

    This function saves the trained RandomForestClassifier model, LabelEncoder, and feature list as a dictionary
    in a pickle file so they can be loaded together later for making predictions.

    Args:
        model (RandomForestClassifier): The trained RandomForestClassifier model.
        label_encoder (LabelEncoder): The LabelEncoder used for encoding the target variable (i.e. alg. config.).
        features (list): The list of feature column names.
        output_file (str): Path to save the model, encoders, and features.
    """
    combined_data = {
        'model': model,
        'label_encoder': label_encoder,
        'features': features
    }
    print(combined_data)
    with open(output_file, 'wb') as f:
        pickle.dump(combined_data, f)

def main():
    """
    Main execution function.

    This function loads the live info and tuning result data from directories, preprocesses it, trains the model,
    and saves both the model and label encoder in a single file. It prints completion messages after training and saving.
    """
    parser = argparse.ArgumentParser(description="Train a RandomForestClassifier on merged AutoPas data.")
    parser.add_argument('--results-dir', type=str, required=True, help="Directory containing live info and tuning"
                                                                       " results CSV files.")
    parser.add_argument('--features', type=str, nargs='+',
                        default=['meanParticlesPerCell', 'maxParticlesPerCell', 'relativeParticlesPerCellStdDev',
                                 'threadCount', 'relativeParticlesPerBlurredCellStdDev', 'skin'],
                        help="List of feature columns for training.")
    parser.add_argument('--test-size', type=float, default=0.2, help="Test set size as a fraction (default: 0.2).")
    parser.add_argument('--n-estimators', type=int, default=100,
                        help="Number of estimators (trees) in the RandomForest (default: 100).")
    parser.add_argument('--output-file', type=str, default='model.pkl', help="Output file to save models and encoders.")
    args = parser.parse_args()

    print("Loading and merging data...")
    merged_df = load_data_from_directory(args.live_info_dir, args.results_dir)

    print("Preprocessing data...")
    X, y, label_encoder = preprocess_data(merged_df, args.features)

    print("Training model...")
    model = train_model(X, y, args.test_size, args.n_estimators)

    print(f"Saving model and encoders to {args.output_file}...")
    save_model_and_encoder(model, label_encoder, args.features, args.output_file)

    print("Model training and saving complete.")


if __name__ == "__main__":
    main()
