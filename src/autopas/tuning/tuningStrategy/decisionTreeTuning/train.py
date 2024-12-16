# @File: train.py
# @Author Abdulkadir Pazar
# @Date 10-08-2024

import os
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from sklearn.multioutput import MultiOutputClassifier
import pickle
import numpy as np

def extract_identifier(filename: str, prefix: str) -> str:
    """
    Extract the timestamp from the filename, assuming it is the last part of the filename before the extension and
    is seperated by underscores as is the naming convention for logging in Autopas.

    Args:
        filename (str): The filename from which to extract the identifier.
        prefix (str): The prefix to locate before the identifier.

    Returns:
        str: The extracted identifier or None if the prefix is not found.
    """
    if filename.startswith(prefix):
        return filename[len(prefix):].split('.')[0]
    return None

def load_data_from_directory(live_info_dir: str, tuning_results_dir: str) -> pd.DataFrame:
    """
    Load and merge multiple live info and tuning results CSV files based on matching timestamps in filenames.

    This function matches files from the live info and tuning results directories by their identifiers and merges
    the data based on the 'Iteration' column.

    Args:
        live_info_dir (str): The directory containing live info CSV files.
        tuning_results_dir (str): The directory containing tuning results CSV files.

    Returns:
        pd.DataFrame: A merged DataFrame containing the combined live info and tuning results from all files.
    """
    merged_dfs = []
    live_info_prefix = "AutoPas_liveInfoLogger_"
    tuning_results_prefix = "AutoPas_tuningResults_"
    live_info_files = os.listdir(live_info_dir)
    tuning_results_files = os.listdir(tuning_results_dir)
    # Create dictionaries of identifiers and filenames for live info and tuning results
    live_info_dict = {
        extract_identifier(f, live_info_prefix): f
        for f in live_info_files if f.startswith(live_info_prefix)
    }
    tuning_results_dict = {
        extract_identifier(f, tuning_results_prefix): f
        for f in tuning_results_files if f.startswith(tuning_results_prefix)
    }
    # Find common identifiers and merge the data
    common_identifiers = set(live_info_dict.keys()) & set(tuning_results_dict.keys())
    for identifier in common_identifiers:
        live_info_path = os.path.join(live_info_dir, live_info_dict[identifier])
        tuning_results_path = os.path.join(tuning_results_dir, tuning_results_dict[identifier])
        live_info_df = pd.read_csv(live_info_path)
        tuning_results_df = pd.read_csv(tuning_results_path)

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

def save_model_and_encoders(model: MultiOutputClassifier, label_encoders: dict, output_file: str) -> None:
    """
    Save the trained model and LabelEncoders to a file.

    This function saves the trained MultiOutputClassifier model and LabelEncoders as a dictionary in a pickle file
    so they can be loaded together later for making predictions.

    Args:
        model (MultiOutputClassifier): The trained MultiOutputClassifier model.
        label_encoders (dict): The LabelEncoders used for encoding target variables.
        output_file (str): Path to save the model and encoders.
    """
    combined_data = {
        'model': model,
        'label_encoders': label_encoders
    }
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
    parser.add_argument('--tuning-results-dir', type=str, required=True, help="Directory containing tuning results CSV files.")
    parser.add_argument('--features', type=str, nargs='+', default=['avgParticlesPerCell', 'maxParticlesPerCell', 'homogeneity', 'maxDensity', 'particlesPerCellStdDev', 'threadCount'], help="List of feature columns for training.")
    parser.add_argument('--targets', type=str, nargs='+', default=['Container', 'Traversal', 'Data Layout', 'Newton 3'], help="List of target columns to predict.")
    parser.add_argument('--test-size', type=float, default=0.2, help="Test set size as a fraction (default: 0.2).")
    parser.add_argument('--n-iterations', type=int, default=100, help="Number of iterations for RandomForest (default: 100).")
    parser.add_argument('--output-file', type=str, default='model.pkl', help="Output file to save models and encoders.")
    args = parser.parse_args()

    print("Loading and merging data...")
    merged_df = load_data_from_directory(args.live_info_dir, args.tuning_results_dir)

    print("Preprocessing data...")
    X, y, label_encoders = preprocess_data(merged_df, args.features, args.targets)

    print("Training model...")
    model = train_model(X, y, args.test_size, args.n_iterations)

    print(f"Saving model and encoders to {args.output_file}...")
    save_model_and_encoders(model, label_encoders, args.output_file)

    print("Model training and saving complete.")

if __name__ == "__main__":
    main()
