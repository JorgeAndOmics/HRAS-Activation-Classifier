import numpy as np
import os
import warnings
import pickle
from itertools import combinations
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import pandas as pd
import torch



class Datahub:

    @staticmethod
    def directory_generator():
        if not os.path.exists(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors'))):
            os.mkdir(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors')))
        if not os.path.exists(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'models'))):
            os.mkdir(os.path.join(os.path.join(os.path.dirname(__file__), 'temp_dir', 'models')))

    @staticmethod
    def warning_silencer():
        warnings.filterwarnings('ignore', message='.*invalid value.*')


    alignment_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'alignment_dataframe.csv'), index_col=0)
    features_df = {
        'w': pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'wcn_dataframe.csv'), index_col=0),
        'p': pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'propka_dataframe.csv'), index_col=0),
        't': pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'pyrosetta_total_energy_dataframe.csv'), index_col=0),
        's':pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'pyrosetta_solvation_dataframe.csv'), index_col=0),
        'e':pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'pyrosetta_electrostatics_dataframe.csv'), index_col=0),
        'v1': pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'pyrosetta_vanderwaals1_dataframe.csv'), index_col=0),
        'v2': pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes' , 'pyrosetta_vanderwaals2_dataframe.csv'), index_col=0)
                   }
    features_df_normalized = None
    features_tensor_normalized = None
    structure_identifier = None
    train_losses = []
    train_accuracies = []
    test_losses = []
    test_accuracies = []
    confusion_matrices = []
    roc_curves = []
    auc_scores = []


def normalizer():
    # Extracting the protein structure codes for later use
    protein_codes = Datahub.alignment_df.index.to_list()
    print(len(protein_codes))

    # Z-score normalization
    scaler = StandardScaler()
    features_normalized = {feature: (scaler.fit_transform(df), df.columns) for feature, df in Datahub.features_df.items()}

    # Impute NaN values with 0 after normalization
    imputer = SimpleImputer(missing_values=np.nan, strategy='constant', fill_value=0)
    features_normalized_imputed = {feature: (imputer.fit_transform(df), columns) for feature, (df, columns) in features_normalized.items()}

    # Convert back to DataFrame for easy NaN handling
    features_normalized_df = {feature: pd.DataFrame(df, index=protein_codes, columns=columns) for feature, (df, columns) in features_normalized_imputed.items()}
    Datahub.features_df_normalized = features_normalized_df

    # Generate csv files for each feature with a comprehension
    for feature, df in features_normalized_df.items():
        df.to_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', f'{feature}_dataframe_normalized.csv'))
        print('CSV file successfully generated for feature at temporary directory')

    return features_normalized_df


def tensor_generator():
    def structure_identifier_generator():
        # Generate a dictionary of positions to structure ids in the ordered alignment dataframe
        structure_identifier = {}
        for i, structure_id in enumerate(Datahub.alignment_df.sort_index().index):
            structure_identifier[i] = structure_id
        Datahub.structure_identifier = structure_identifier
        with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'structure_identifier.pickle'), "wb") as f:
            pickle.dump(structure_identifier, f)
            print('Structure identifier successfully dumped to pickle file')

        return structure_identifier

    # Create all combinations of normalized tensors
    for r in range(1, len(Datahub.features_df_normalized) + 1):
        for subset in combinations(Datahub.features_df_normalized.items(), r):
            # Sort each dataframe by index, convert to numpy array, and stack to create a tensor
            feature_arrays_normalized = [df.sort_index().values for name, df in subset]
            feature_tensor_normalized = np.dstack(feature_arrays_normalized)
            feature_tensor_normalized = torch.from_numpy(feature_tensor_normalized).float()
            print(feature_tensor_normalized.shape)

            # Create a name for each combination for saving purposes
            subset_names = [name for name, df in subset]
            combination_name = ''.join(subset_names)

            # Save each tensor to a pickle file
            if feature_tensor_normalized.shape[2] != 7:
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors', f'{combination_name}_tensor.pickle'), "wb") as f:
                    pickle.dump(feature_tensor_normalized, f)
                    print(f'{combination_name} normalized features tensor successfully dumped to pickle file')
            else:
                with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors', 'complete_tensor.pickle'), "wb") as f:
                    pickle.dump(feature_tensor_normalized, f)
                    print('Complete normalized features tensor successfully dumped to pickle file')

    structure_identifier_generator()

    return feature_tensor_normalized if feature_tensor_normalized.shape[2] == 7 else None




if __name__ == '__main__':

    Datahub.warning_silencer()
    Datahub.directory_generator()

    normalizer = normalizer()

    complete_tensor = tensor_generator()

    # Load the data here
    # with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'tensors', 'complete_tensor.pickle'), 'rb') as f:
    #     tensor_data = pickle.load(f)
    # metadata_df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'temp_dir', 'result_dataframes', 'metadata_dataframe.csv'), index_col=0)
    # with open(os.path.join(os.path.dirname(__file__), 'temp_dir', 'pickles', 'structure_identifier.pickle'), "rb") as f:
    #     structure_identifier = pickle.load(f)


