# @File: predict.py
# @Author Abdulkadir Pazar
# @Date 10-08-2024

import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
import pickle


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


def load_data_from_directory(live_info_dir: str, tuning_results_dir: str) -> pd.DataFrame:
    """
    Load and merge multiple live info and tuning results CSV files based on matching timestamps in filenames.

    This function matches files from the live info and tuning results directories by their timestamps and merges
    the data based on the 'Iteration' column.

    Args:
        live_info_dir (str): The directory containing live info CSV files.
        tuning_results_dir (str): The directory containing tuning results CSV files.

    Returns:
        pd.DataFrame: A merged DataFrame containing the combined live info and tuning results from all files.
    """
    merged_dfs = []

    # Get all the filenames in both directories
    live_info_files = os.listdir(live_info_dir)
    tuning_results_files = os.listdir(tuning_results_dir)

    # Create a dictionary to map timestamps to file paths for each directory
    live_info_dict = {extract_timestamp(f): f for f in live_info_files}
    tuning_results_dict = {extract_timestamp(f): f for f in tuning_results_files}

    # Find common timestamps between live info and tuning results files
    common_timestamps = set(live_info_dict.keys()) & set(tuning_results_dict.keys())

    # Iterate over common timestamps and load/merge corresponding files
    for timestamp in common_timestamps:
        live_info_path = os.path.join(live_info_dir, live_info_dict[timestamp])
        tuning_results_path = os.path.join(tuning_results_dir, tuning_results_dict[timestamp])

        # Load the live info and tuning results
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


def preprocess_data(merged_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Preprocess the merged DataFrame for model training.

    This function selects the relevant features from the live info and the target variables from the tuning results.
    It encodes categorical target variables using LabelEncoder and splits the merged DataFrame into features (X)
    and target labels (y).

    Args:
        merged_df (pd.DataFrame): The merged DataFrame containing live info and tuning results.

    Returns:
        X (pd.DataFrame): A DataFrame of feature columns used for training.
        y (pd.DataFrame): A DataFrame of target columns to be predicted.
        label_encoders (dict): A dictionary of LabelEncoders used for encoding the target variables.
    """
    # Select the features for training
    features = ['avgParticlesPerCell', 'maxParticlesPerCell', 'homogeneity', 'maxDensity', 'particlesPerCellStdDev',
                'threadCount']

    # Select the target
    targets = ['Container', 'Traversal', 'Load Estimator', 'Data Layout', 'Newton 3', 'CellSizeFactor']

    # Encode using LabelEncoder
    label_encoders = {target: LabelEncoder() for target in targets}

    for target in targets:
        merged_df[target] = label_encoders[target].fit_transform(merged_df[target])

    X = merged_df[features]
    y = merged_df[targets]

    return X, y, label_encoders


def train_model(X: pd.DataFrame, y: pd.DataFrame) -> dict[str, RandomForestClassifier]:
    """
    Train a RandomForestClassifier for each target variable.

    This function splits the data into training and testing sets and trains a RandomForestClassifier for each target
    variable (e.g., 'Container', 'Traversal', etc.). It also evaluates each model by calculating the accuracy on
    the test data.

    Args:
        X (pd.DataFrame): The feature set for training the model.
        y (pd.DataFrame): The target set for training the model.

    Returns:
        models (dict): A dictionary of trained RandomForestClassifier models for each target variable.
    """
    # Split data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Train a RandomForestClassifier for each target
    models = {}
    for i, column in enumerate(y.columns):
        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train, y_train.iloc[:, i])
        models[column] = model

        # Test accuracy on test set for this model
        y_pred = model.predict(X_test)
        accuracy = accuracy_score(y_test.iloc[:, i], y_pred)
        print(f'Accuracy for {column}: {accuracy:.2f}')

    return models


def save_model_and_encoders(models: dict[str, RandomForestClassifier], label_encoders: dict[str, LabelEncoder]) -> None:
    """
    Save the trained models and LabelEncoders to a single file.

    This function saves the trained models and LabelEncoders as a dictionary in a single pickle file
    so they can be loaded together later for making predictions.

    Args:
        models (dict): The trained models to be saved.
        label_encoders (dict): The label encoders used for encoding target variables.
    """
    # Combine models and label_encoders in a single dictionary
    combined_data = {
        'models': models,
        'label_encoders': label_encoders
    }

    # Save combined data to a single pickle file
    with open('model.pkl', 'wb') as f:
        pickle.dump(combined_data, f)


# Main function
if __name__ == "__main__":
    """
    Main execution function.

    This function loads the live info and tuning result data from directories, preprocesses it, trains the models,
    and saves both the models and label encoders in a single file. It prints completion messages after training and saving.
    """
    # Directories for live info and tuning results
    live_info_dir = 'data/live_info'
    tuning_results_dir = 'data/tuning_results'

    # Load and preprocess the data
    merged_df = load_data_from_directory(live_info_dir, tuning_results_dir)
    X, y, label_encoders = preprocess_data(merged_df)

    # Train the model
    models = train_model(X, y)

    # Save the trained models and encoders in a file
    save_model_and_encoders(models, label_encoders)

    print("Model training and saving complete.")

