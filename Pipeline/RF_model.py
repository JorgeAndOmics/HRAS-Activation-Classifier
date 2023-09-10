import pickle
import os
import numpy as np
import pandas as pd
import argparse
from sklearn.manifold import TSNE
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import (accuracy_score, classification_report, confusion_matrix,
                             roc_curve, roc_auc_score, precision_recall_curve, log_loss,
                             cohen_kappa_score, matthews_corrcoef, brier_score_loss)


def directory_generator():
    if not os.path.exists(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'models', 'RF'))):
        os.mkdir(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'models', 'RF')))

def load_tensors(directory):
    """
    Load all tensor files from the specified directory.

    Args:
    - directory (str): Path to the directory containing tensor files.

    Returns:
    - List of tensors and corresponding filenames.
    """
    tensor_files = [f for f in os.listdir(directory) if f.endswith('.pickle')]
    tensors = []
    filenames = []

    for tensor_file in tensor_files:
        with open(os.path.join(directory, tensor_file), 'rb') as file:
            tensor = pickle.load(file)
            tensors.append(tensor)
            filenames.append(tensor_file.split(".")[0])

    print(f"Loaded {len(tensors)} tensors.")
    return tensors, filenames


def apply_tsne(tensors, filenames, seed):
    """
    Apply t-SNE to tensors to reduce dimensionality to 2D.

    Args:
    - tensors (list): List of tensors.
    - filenames (list): List of tensor filenames.

    Returns:
    - 2D numpy array after t-SNE transformation.
    """
    tsne = TSNE(n_components=2, random_state=seed)
    transformed_data = []

    for tensor, filename in zip(tensors, filenames):
        print(f'Applying t-SNE transformation to {filename}...')
        flat_data = tensor.reshape(tensor.shape[0], -1)
        embedded = tsne.fit_transform(flat_data)
        transformed_data.append(embedded)
        print('t-SNE transformation completed.')

    return np.vstack(transformed_data)

def train_random_forest(X, y, seed):
    """
    Train a random forest classifier on the data and compute various metrics.

    Args:
    - X (numpy array): Input data.
    - y (list): Numeric labels (1 for 'active', 0 for 'inactive').

    Returns:
    - Trained classifier, metrics dictionary containing accuracy, classification report,
      confusion matrix, and additional metrics.
    """

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=seed)

    # Define the parameter grid
    param_grid = {
        'n_estimators': [10, 50, 100, 200, 500],
        'max_depth': [None, 10, 20, 30, 40, 50],
        'min_samples_split': [2, 5, 10, 20],
        'min_samples_leaf': [1, 2, 4, 8, 16],
        'max_features': ['sqrt', 'log2', None],
        'criterion': ['gini', 'entropy']
    }

    # Create the base model to tune
    rf = RandomForestClassifier(random_state=seed)

    # Instantiate the grid search model
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid,
                               cv=3, n_jobs=-1, verbose=2, scoring='accuracy')

    # Fit the grid search model
    grid_search.fit(X_train, y_train)

    # Get the best estimator
    best_clf = grid_search.best_estimator_
    best_params = grid_search.best_params_

    y_pred = best_clf.predict(X_test)
    y_pred_prob = best_clf.predict_proba(X_test)[:, 1]

    metrics = {
        'accuracy': accuracy_score(y_test, y_pred),
        'classification_report': classification_report(y_test, y_pred),
        'confusion_matrix': confusion_matrix(y_test, y_pred),
        'roc_curve': roc_curve(y_test, y_pred_prob),
        'roc_auc': roc_auc_score(y_test, y_pred_prob),
        'precision_recall_curve': precision_recall_curve(y_test, y_pred_prob),
        'log_loss': log_loss(y_test, y_pred_prob),
        'cohen_kappa': cohen_kappa_score(y_test, y_pred),
        'matthews_corrcoef': matthews_corrcoef(y_test, y_pred),
        'brier_score': brier_score_loss(y_test, y_pred_prob),
        'feature_importances': best_clf.feature_importances_,
        'best_parameters': best_params  # Storing the best parameters
    }

    print(f'Training completed with accuracy: {metrics["accuracy"]:.4f}')

    return best_clf, metrics


def extract_and_store_results(classifier, metrics, output_dir, filename):
    """
    Store results and metrics to pickle files.

    Args:
    - classifier (RandomForestClassifier): Trained random forest classifier.
    - metrics (dict): Dictionary containing various metrics.
    - output_dir (str): Directory to store results.
    """
    with open(os.path.join(output_dir, f'RF-{filename.split("_")[0]}_model.pkl'), 'wb') as f:
        pickle.dump(classifier, f)

    with open(os.path.join(output_dir, f'RF-{filename.split("_")[0]}_model_metrics.pkl'), 'wb') as f:
        pickle.dump(metrics, f)  # Store the entire metrics dictionary

    print('Results stored successfully.')

def main(tensor, labels, filename, seed, index):
    """
    Main execution for each tensor.

    Args:
    - tensor (numpy array): The tensor data.
    - labels (list): Activity status labels.
    - filename (str): Tensor filename (used for logging).
    - seed (int): Random seed for reproducibility.

    Returns:
    - None
    """

    print(f'Processing {filename}...[{index}/{len(filenames)}]')

    # Apply t-SNE
    transformed_data = apply_tsne([tensor], [filename], seed)

    # Train random forest and get metrics
    print(f'Training random forest classifier on {filename}...')
    clf, metrics = train_random_forest(transformed_data, labels, seed)

    # Extract and store results
    result_dir = os.path.join(os.path.dirname(__file__), 'temp_dir', 'models', 'RF', filename)
    os.makedirs(result_dir, exist_ok=True)
    extract_and_store_results(clf, metrics, result_dir, filename)

    print(f"Processing for {filename} completed.\n")



if __name__ == "__main__":

    # SCRIPT CALL
    parser = argparse.ArgumentParser(description="Script for analysis configuration")
    parser.add_argument("--seed", type=int, default=1,
                        help="Seed for analysis (default: 1)")
    args = parser.parse_args()
    seed = args.seed

    # Generate directory
    directory_generator()

    # RANDOM FOREST
    # Load tensors
    print('Initializing Random Forest model generation...')
    tensors, filenames = load_tensors(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors'))

    # Load and convert labels
    metadata_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'metadata_dataframe.csv'))
    labels = metadata_df["Read Activity Status"].replace({'inactive': 0, 'active': 1}).tolist()

    for index, (tensor, filename) in enumerate(zip(tensors, filenames), start=1):
        main(tensor, labels, filename, seed, index)
