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
                             model_prefix: str, save_python: bool, save_treelite: bool) -> None:
    """
    Save training artifacts for Python and/or Treelite-based prediction.

    When Python output is enabled, this function saves the trained RandomForestClassifier models and LabelEncoders
    for both pairwise and triwise interactions, together with the list of features, into a single pickle bundle.
    This bundle is loaded later for making predictions.

    When Treelite output is enabled, this function exports Treelite checkpoint files, classes files, and one
    features file. These files are later loaded separately: the model checkpoint for inference,
    the classes file for decoding class ids to labels, and the features file for consistent input order.

    Args:
        pairwise_model (RandomForestClassifier): The trained RandomForestClassifier model for pairwise interactions.
        triwise_model (RandomForestClassifier): The trained RandomForestClassifier model for triwise interactions.
        pairwise_label_encoder (LabelEncoder): The LabelEncoder used for encoding the target variable (i.e. alg.
        config.).
        triwise_label_encoder (LabelEncoder): The LabelEncoder used for encoding the target variable (i.e. alg.
        config.).
        features (list): The list of feature column names.
        model_prefix (str): Base name/prefix used for generated output files.
        save_python (bool): If True, writes '<model-prefix>.pkl'.
        save_treelite (bool): If True, writes Treelite outputs:
            - '<model-prefix>_pairwise.tl' and '<model-prefix>_pairwise_classes.txt' (if pairwise model exists),
            - '<model-prefix>_triwise.tl' and '<model-prefix>_triwise_classes.txt' (if triwise model exists),
            - '<model-prefix>_features.json' (if at least one Treelite model exists).
    """
    def _treelite_skip_reason(model, encoder) -> str | None:
        """
        Return the reason why Treelite export should be skipped for a model.

        Treelite export is skipped if no model/encoder exists for an interaction type, or if the trained
        classifier has less than two classes.

        Args:
            model: Trained classifier model for one interaction type.
            encoder: LabelEncoder fitted for the same interaction type.

        Returns:
            str | None: Human-readable skip reason, or None if export may proceed.
        """
        if model is None or encoder is None:
            return "no tuning results found"
        
        # Treelite import for sklearn classifiers requires at least two classes.
        # A single-class model can still be used in the Python/pickle path, but it will always
        # predict the same label and is therefore not useful for tuning decisions.
        # Such models cannot be exported to Treelite (.tl).
        if getattr(model, "n_classes_", 0) < 2:
            return ("trained model has < 2 classes; Treelite export requires at least 2 classes. "
                    "Collect more training data for this interaction type.")

    # Save to .pkl
    if save_python:
        output_file = f"{model_prefix}.pkl"
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
        print(f"Saved pickle bundle: {output_file}")

    tl_created = False

    if save_treelite:
        # Save Treelite pairwise model and classes using shared model prefix.
        pairwise_output_file = f"{model_prefix}_pairwise.tl"
        pairwise_classes_file = f"{model_prefix}_pairwise_classes.txt"

        reason = _treelite_skip_reason(pairwise_model, pairwise_label_encoder)
        if reason is not None:
            print(f"Skipping {pairwise_output_file}: {reason}")
        else:
            tl_model = treelite.sklearn.import_model(pairwise_model)
            tl_model.serialize(pairwise_output_file)
            np.savetxt(pairwise_classes_file, pairwise_label_encoder.classes_, fmt="%s")
            print(f"Saved Treelite model: {pairwise_output_file}")
            print(f"Saved classes: {pairwise_classes_file}")
            tl_created = True

        # Save Treelite triwise model and classes using shared model prefix.
        triwise_output_file = f"{model_prefix}_triwise.tl"
        triwise_classes_file = f"{model_prefix}_triwise_classes.txt"

        reason = _treelite_skip_reason(triwise_model, triwise_label_encoder)
        if reason is not None:
            print(f"Skipping {triwise_output_file}: {reason}")
        else:
            tl_model = treelite.sklearn.import_model(triwise_model)
            tl_model.serialize(triwise_output_file)
            np.savetxt(triwise_classes_file, triwise_label_encoder.classes_, fmt="%s")
            print(f"Saved Treelite model: {triwise_output_file}")
            print(f"Saved classes: {triwise_classes_file}")
            tl_created = True

    # Dump <model-prefix>_features.json only if at least one treelite model was created
    if tl_created:
        features_file = f"{model_prefix}_features.json"
        with open(features_file, "w") as f:
            json.dump(list(features), f)
        print(f"Saved features: {features_file}")


def main():
    """
    Main execution function.

    This function loads the live info and tuning result data from directories, preprocesses it, trains the model,
    and saves the trained outputs in Python pickle and/or Treelite format depending on the selected options.
    It prints completion messages after training and saving.
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
    parser.add_argument('--model-prefix', type=str, default='model',
                        help="Base name/prefix for outputs. Produces '<model-prefix>.pkl' (Python) and Treelite related files like '<model-prefix>_pairwise.tl'.")
    parser.add_argument('--save-python', action=argparse.BooleanOptionalAction, default=True,
                        help="Enable or disable writing '<model-prefix>.pkl' via '--save-python' / '--no-save-python' "
                             "(default: enabled).")
    parser.add_argument('--save-treelite', action=argparse.BooleanOptionalAction, default=True,
                        help="Enable or disable Treelite outputs via '--save-treelite' / '--no-save-treelite' "
                             "using '<model-prefix>' (default: enabled).")
    args = parser.parse_args()

    # Use one prefix for all output files.
    model_prefix = args.model_prefix
    if not model_prefix:
        parser.error("--model-prefix must not be empty.")

    # Check that at least one output target is enabled.
    if not args.save_python and not args.save_treelite:
        parser.error("At least one output target must be enabled: --save-python and/or --save-treelite.")

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

    print(f"Saving model outputs with prefix '{model_prefix}' "
          f"(save-python={args.save_python}, save-treelite={args.save_treelite})...")
    save_models_and_encoders(model_pairwise, model_triwise, label_encoder_pairwise, label_encoder_triwise,
                             args.features, model_prefix, args.save_python, args.save_treelite)

    print("Model training and saving complete.")


if __name__ == "__main__":
    main()
