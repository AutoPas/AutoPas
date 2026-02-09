# @File: train.py
# @Author Abdulkadir Pazar
# @Date 10-08-2024

import os
import argparse
import glob

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from sklearn.multioutput import MultiOutputClassifier
import pickle
import numpy as np

import json
import treelite
import treelite.sklearn


def extract_filename(path_to_file: str) -> str:
    """
    Extract the filename from a path.

    E.g. './path/to/file.csv' becomes 'file.csv'.

    Args:
        path_to_file (str): The path to the file.
    Returns:
         str: The filename.
    """
    return os.path.basename(path_to_file)

def extract_identifier(filename: str) -> str:
    """
    Extract the rank and timestamp from the filename, to serve as a unique identifier to match Live Info csvs to
    tuning results csvs.

    E.g. 'AutoPas_liveInfoLogger_Rank0_2025-05-13_14-11-48.csv' becomes 'Rank0_2025-05-13_14-11-48.csv'
    'AutoPas_tuningResults_Rank0_pairwise_2025-05-13_14-11-48.csv' becomes 'Rank0_pairwise_2025-05-13_14-11-48.csv'
    (and therefore the 'pairwise' part must then be removed by remove_interaction_type so they match).

    Args:
        filename (str): The filename from which to extract the identifier.

    Returns:
        str: The extracted identifier.
    """
    return filename.split('_', 2)[2]

def remove_interaction_type(identifier: str) -> str:
    """
    Remove the interaction type from the identifier produced by extract_identifier applied to a TuningResults csv.

    E.g. 'Rank0_pairwise_2025-05-13_14-11-48.csv' becomes 'Rank0_2025-05-13_14-11-48.csv'

    Args:
        identifier (str): The filename of a tuning results file, after extract_identifier has been applied to it.

    Returns:
        str: The identifier with the interaction type removed.
    """
    parts = identifier.split('_')
    return '_'.join(parts[:1] + parts[2:])


def load_data_from_directory(results_dir: str) -> tuple:
    """
    Load and merge multiple live info and tuning results CSV files assuming leaf folders contain all experiments with
    each folder containing matching live info and tuning result CSV files. Files are matched together using their rank
    and timestamps.

    Args:
        results_dir (str): The directory containing all experiment CSV files.

    Returns:
        final_merged_df_pairwise (pd.DataFrame): A merged DataFrame containing the combined live info and pairwise tuning results from all files.
        final_merged_df_triwise (pd.DataFrame): A merged DataFrame containing the combined live info and triwise tuning results from all files.
    """
    merged_dfs_pairwise = []
    merged_dfs_triwise = []

    # Loop through all leaf folders i.e. the ones in which the csv files are outputted in

    for root, dirs, files in os.walk(results_dir):
        if not dirs:
            # Get Live Info files and sort them
            live_info_file_pattern = os.path.join(root,"*liveInfoLogger*")
            live_info_file_list = glob.glob(live_info_file_pattern)

            # Get Tuning Results files
            tuning_results_file_pattern = os.path.join(root,"*tuningResults*")
            tuning_results_file_list = glob.glob(tuning_results_file_pattern)

            # Separate Tuning Results into Pairwise and Triwise
            tuning_results_file_list_pairwise = [f for f in tuning_results_file_list if '_pairwise' in f]
            tuning_results_file_list_triwise = [f for f in tuning_results_file_list if '_triwise' in f]

            # Create dictionaries mapping identifiers to files and use it to get a list of matching file pairs
            live_info_map = {extract_identifier(extract_filename(f)): f for f in live_info_file_list}
            tuning_results_pairwise_map = {remove_interaction_type(extract_identifier(extract_filename(f))):
                                               f for f in tuning_results_file_list_pairwise}
            tuning_results_triwise_map = {remove_interaction_type(extract_identifier(extract_filename(f))):
                                              f for f in tuning_results_file_list_triwise}

            common_identifiers_pairwise = set(live_info_map) & set(tuning_results_pairwise_map)
            common_identifiers_triwise = set (live_info_map) & set(tuning_results_triwise_map)

            csv_file_pairs_pairwise = [(live_info_map[id], tuning_results_pairwise_map[id]) for id in common_identifiers_pairwise]
            csv_file_pairs_triwise = [(live_info_map[id], tuning_results_triwise_map[id]) for id in common_identifiers_triwise]

            # Pairwise
            for [live_info_file, tuning_results_file] in csv_file_pairs_pairwise:
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
                merged_dfs_pairwise.append(merged_df)

            # Triwise
            for [live_info_file, tuning_results_file] in csv_file_pairs_triwise:
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
                merged_dfs_triwise.append(merged_df)

    # Combine all merged DataFrames into one
    if len(merged_dfs_pairwise) > 0:
        print(f"Merging {len(merged_dfs_pairwise)} CSV files for the pairwise training data")
        final_merged_df_pairwise = pd.concat(merged_dfs_pairwise, ignore_index=True)
    else:
        final_merged_df_pairwise = None
    if len(merged_dfs_triwise) > 0:
        final_merged_df_triwise = pd.concat(merged_dfs_triwise, ignore_index=True)
    else:
        final_merged_df_triwise = None

    return final_merged_df_pairwise, final_merged_df_triwise


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

    merged_df['alg_config'] = label_encoder.fit_transform(merged_df['alg_config'])
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
    model = RandomForestClassifier(n_estimators=n_estimators, random_state=42, max_features=1.)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    accuracy = accuracy_score(y_test, y_pred)
    print(f'Accuracy: {accuracy:.3f}')

    return model


def save_models_and_encoders(pairwise_model: RandomForestClassifier, triwise_model: RandomForestClassifier,
                             pairwise_label_encoder: dict, triwise_label_encoder: dict, features: list,
                             output_files: list) -> None:
    """
    Save the trained models, LabelEncoder, and feature list to a file.

    This function saves the trained RandomForestClassifier models and LabelEncoders for both pairwise and triwise
    interactions, as well as the list of features used, as a dictionary in a pickle file so they can be loaded together
    later for making predictions.

    Additionally, it attempts to export Treelite (.tl) models for pairwise and/or triwise interactions. For each successful
    export, a matching '<model>_classes.txt' file with label encoder's class mapping is created. If at least one Treelite 
    model is successfully created, the feature list is written once to `features.json`.

    Args:
        pairwise_model (RandomForestClassifier): The trained RandomForestClassifier model for pairwise interactions.
        triwise_model (RandomForestClassifier): The trained RandomForestClassifier model for triwise interactions.
        pairwise_label_encoder (LabelEncoder): The LabelEncoder used for encoding the target variable (i.e. alg.
        config.).
        triwise_label_encoder (LabelEncoder): The LabelEncoder used for encoding the target variable (i.e. alg.
        config.).
        features (list): The list of feature column names.
        output_file (list): The list of paths to save the models, encoders, and features.
    """

    def _is_ext(name: str, ext: str) -> bool:
        return name.lower().endswith(ext)

    def _strip_ext(name: str, ext: str) -> str:
        return name[:-len(ext)] if _is_ext(name, ext) else name

    def _treelite_skip_reason(model, encoder) -> str | None:
        if model is None or encoder is None:
            return "no tuning results found"
        if getattr(model, "n_classes_", 0) < 2:
            return "trained model has < 2 classes"

    # Save to .pkl
    for output_file in output_files:
        if not _is_ext(output_file, ".pkl"):
            continue
        
        combined_data = {
            'pairwise_model': pairwise_model,
            'triwise_model': triwise_model,
            'pairwise_label_encoder': pairwise_label_encoder,
            'triwise_label_encoder': triwise_label_encoder,
            'features': features
        }
        print(combined_data)
        with open(output_file, 'wb') as f:
            pickle.dump(combined_data, f)

    tl_created = False

    # Save to .tl and .txt
    for output_file in output_files:
        if not _is_ext(output_file, ".tl"):
            continue

        base = _strip_ext(output_file, ".tl")
        classes_file = f"{base}_classes.txt"

        if "pairwise" in base.lower():
            reason = _treelite_skip_reason(pairwise_model, pairwise_label_encoder)
            if reason is not None:
                print(f"Skipping {output_file}: {reason}")
                continue

            tl_model = treelite.sklearn.import_model(pairwise_model)
            tl_model.serialize(output_file)

            np.savetxt(classes_file, pairwise_label_encoder.classes_, fmt="%s")
            
            print(f"Saved model: {output_file}")
            print(f"Saved classes: {classes_file}")

            tl_created = True

        elif "triwise" in base.lower():
            reason = _treelite_skip_reason(triwise_model, triwise_label_encoder)
            if reason is not None:
                print(f"Skipping {output_file}: {reason}")
                continue

            tl_model = treelite.sklearn.import_model(triwise_model)
            tl_model.serialize(output_file)

            np.savetxt(classes_file, triwise_label_encoder.classes_, fmt="%s")
            
            print(f"Saved model: {output_file}")
            print(f"Saved classes: {classes_file}")

            tl_created = True

        else:
            print(f"Skipping .tl output '{output_file}': filename must contain 'pairwise' or 'triwise'.")

    # Dump features.json only if at least one treelite model was created
    if tl_created:
        with open("features.json", "w") as f:
            json.dump(list(features), f)
        print("Saved features: features.json")


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
                        default=['meanParticlesPerCell', 'medianParticlesPerCell', 'maxParticlesPerCell', 'relativeParticlesPerCellStdDev',
                                 'threadCount', 'numCells', 'numEmptyCells', 'skin'],
                        help="List of feature columns for training.")
    parser.add_argument('--test-size', type=float, default=0.1, help="Test set size as a fraction (default: 0.1).")
    parser.add_argument('--n-estimators', type=int, default=100,
                        help="Number of estimators (trees) in the RandomForest (default: 100).")
    parser.add_argument('--output-files', type=str, nargs='+', 
                        default=['model.pkl', 'model_pairwise.tl', 'model_triwise.tl'], 
                        help="Output file to save models and encoders. Supported: .pkl and .tl files.")
    args = parser.parse_args()

    print("Loading and merging data...")
    merged_df_pairwise, merged_df_triwise = load_data_from_directory(args.results_dir)

    if merged_df_pairwise is None:
        print("No tuning results were found for pairwise interactions. Skipping the training of a pairwise model.")
        model_pairwise = None
        label_encoder_pairwise = None
    else:
        print("Preprocessing pairwise data...")
        X, y, label_encoder_pairwise = preprocess_data(merged_df_pairwise, args.features)

        print("Training pairwise model...")
        model_pairwise = train_model(X, y, args.test_size, args.n_estimators)

    if merged_df_triwise is None:
        print("No tuning results were found for triwise interactions. Skipping the training of a triwise model.")
        model_triwise = None
        label_encoder_triwise = None
    else:
        print("Preprocessing triwise data...")
        X, y, label_encoder_triwise = preprocess_data(merged_df_triwise, args.features)

        print("Training triwise model...")
        model_triwise = train_model(X, y, args.test_size, args.n_estimators)


    print(f"Saving models and encoders to {args.output_files}...")
    save_models_and_encoders(model_pairwise, model_triwise, label_encoder_pairwise, label_encoder_triwise, args.features, args.output_files)

    print("Model training and saving complete.")


if __name__ == "__main__":
    main()
