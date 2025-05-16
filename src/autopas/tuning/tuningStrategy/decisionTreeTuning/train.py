# @File: train.py
# @Author Abdulkadir Pazar
# @Date 10-08-2024

import os
import argparse
from pyexpat import features

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from sklearn.multioutput import MultiOutputClassifier
import pickle
import numpy as np

def extract_timestamp(filename: str) -> str:
    """
    Extract the timestamp from the filename, assuming it is the last part of the filename before the extension and
    is seperated by underscores as is the naming convention for logging in Autopas.

    Args:
        filename (str): The filename from which to extract the timestamp.

    Returns:
        str: The extracted timestamp.
    """
    # Split filename by underscores and remove the extension
    return filename.split('_')[-1].split('.')[0]


def load_data_from_directory(results_dir: str) -> pd.DataFrame:
    """
    Load and merge multiple live info and tuning results CSV files assuming leaf folders contain all experiments with
    each folder containing a live info and tuning result CSV file.

    This function matches files from the live info and tuning results directories by their timestamps and merges
    the data based on the 'Iteration' column.

    Args:
        results_dir (str): The directory containing all experiment CSV files.

    Returns:
        pd.DataFrame: A merged DataFrame containing the combined live info and tuning results from all files.
    """
    merged_dfs = []

    # Loop through all leaf folders i.e. the ones in which the csv files are outputted in

    for root, dirs, files in os.walk(results_dir):
        if not dirs:
            # Get Live Info file

            live_info_file_pattern = os.path.join(root,"*liveInfoLogger*")
            live_info_file_list = glob.glob(live_info_file_pattern)

            if len(live_info_file_list) == 0:
                print(f"No Live Info csv file found in folder {root}")

            if len(live_info_file_list) > 1:
                print(f"More than one Live Info csv file found in folder {root}")

            live_info_file = live_info_file_list[0]

            # Get Iteration file

            tuning_results_file_pattern = os.path.join(root,"*tuningResults*")
            tuning_results_file_list = glob.glob(tuning_results_file_pattern)

            if len(tuning_results_file_list) == 0:
                print(f"No Iteration csv file found in folder {root}")

            if len(tuning_results_file_list) > 1:
                print(f"More than one Iteration csv file found in folder {root}")

            tuning_results_file = tuning_results_file_list[0]

            # Load the live info and tuning results
            live_info_df = pd.read_csv(live_info_file)
            tuning_results_df = pd.read_csv(tuning_results_file)

            # Merge them on 'Iteration' column
            merged_df = pd.merge(live_info_df, tuning_results_df, on='Iteration', how='right')

            # Clean column names
            merged_df.columns = merged_df.columns.str.strip()

            # Append the merged DataFrame to the list
            merged_dfs.append(merged_df)

    # Combine all merged DataFrames into one
    final_merged_df = pd.concat(merged_dfs, ignore_index=True)

    return final_merged_df


def preprocess_data(merged_df: pd.DataFrame, features: list, targets: list) -> tuple:
    """
    Preprocess the merged DataFrame for model training.

    This function selects the relevant features from the live info and the target variables from the tuning results.
    It encodes categorical target variables using LabelEncoder and splits the merged DataFrame into features (X)
    and target labels (y).

    Args:
        merged_df (pd.DataFrame): The merged DataFrame containing live info and tuning results.
        features (list): List of column names to be used as features.
        targets (list): List of column names to be used as targets.

    Returns:
        X (pd.DataFrame): A DataFrame of feature columns used for training.
        y (pd.DataFrame): A DataFrame of target columns to be predicted.
        label_encoders (dict): A dictionary of LabelEncoders used for encoding the target variables.
    """
    label_encoders = {target: LabelEncoder() for target in targets}
    for target in targets:
        merged_df[target] = label_encoders[target].fit_transform(merged_df[target])
    X = merged_df[features]
    y = merged_df[targets]
    return X, y, label_encoders


def train_model(X: pd.DataFrame, y: pd.DataFrame, test_size: float, n_iterations: int) -> MultiOutputClassifier:
    """
    Train a MultiOutputClassifier with a RandomForestClassifier for the target vector.

    This function splits the data into training and testing sets and trains a MultiOutputClassifier,
    which supports multi-output classification by wrapping around a RandomForestClassifier.
    It evaluates the model by calculating the accuracy on the test data for each target.

    Args:
        X (pd.DataFrame): The feature set for training the model.
        y (pd.DataFrame): The target set for training the model.
        test_size (float): Fraction of data to be used for testing.
        n_iterations (int): Number of estimators for the RandomForestClassifier.

    Returns:
        MultiOutputClassifier: The trained MultiOutputClassifier model.
    """
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    base_model = RandomForestClassifier(n_estimators=n_iterations, random_state=42)
    model = MultiOutputClassifier(base_model)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    for i, column in enumerate(y.columns):
        accuracy = accuracy_score(y_test.iloc[:, i], y_pred[:, i])
        print(f'Accuracy for {column}: {accuracy:.2f}')

    return model


def save_model_and_encoders(model: MultiOutputClassifier, label_encoders: dict, features: list, output_file: str) -> None:
    """
    Save the trained model, LabelEncoders, and feature list to a file.

    This function saves the trained MultiOutputClassifier model, LabelEncoders, and feature list as a dictionary
    in a pickle file so they can be loaded together later for making predictions.

    Args:
        model (MultiOutputClassifier): The trained MultiOutputClassifier model.
        label_encoders (dict): The LabelEncoders used for encoding target variables.
        features (list): The list of feature column names.
        output_file (str): Path to save the model, encoders, and features.
    """
    combined_data = {
        'model': model,
        'label_encoders': label_encoders,
        'features': features
    }
    print(combined_data)
    with open(output_file, 'wb') as f:
        pickle.dump(combined_data, f)

def main():
    """
    Main execution function.

    This function loads the live info and tuning result data from directories, preprocesses it, trains the model,
    and saves both the model and label encoders in a single file. It prints completion messages after training and saving.
    """
    parser = argparse.ArgumentParser(description="Train a RandomForestClassifier on merged AutoPas data.")
    parser.add_argument('--live-info-dir', type=str, required=True, help="Directory containing live info CSV files.")
    parser.add_argument('--tuning-results-dir', type=str, required=True,
                        help="Directory containing tuning results CSV files.")
    parser.add_argument('--features', type=str, nargs='+',
                        default=['meanParticlesPerCell', 'maxParticlesPerCell', 'relativeParticlesPerCellStdDev',
                                 'threadCount', 'relativeParticlesPerBlurredCellStdDev', 'skin'],
                        help="List of feature columns for training.")
    parser.add_argument('--targets', type=str, nargs='+',
                        default=['Container', 'Traversal', 'Load Estimator', 'Data Layout', 'Newton 3', 'CellSizeFactor'],
                        help="List of target columns to predict.")
    parser.add_argument('--test-size', type=float, default=0.2, help="Test set size as a fraction (default: 0.2).")
    parser.add_argument('--n-iterations', type=int, default=100,
                        help="Number of iterations for RandomForest (default: 100).")
    parser.add_argument('--output-file', type=str, default='model.pkl', help="Output file to save models and encoders.")
    args = parser.parse_args()

    print("Loading and merging data...")
    merged_df = load_data_from_directory(args.live_info_dir, args.tuning_results_dir)

    print("Preprocessing data...")
    X, y, label_encoders = preprocess_data(merged_df, args.features, args.targets)

    print("Training model...")
    model = train_model(X, y, args.test_size, args.n_iterations)

    print(f"Saving model and encoders to {args.output_file}...")
    save_model_and_encoders(model, label_encoders, args.features, args.output_file)

    print("Model training and saving complete.")


if __name__ == "__main__":
    main()
