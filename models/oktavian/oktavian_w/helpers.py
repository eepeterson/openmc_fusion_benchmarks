import numpy as np
import pandas as pd

def process_n_tally(tally_dataframe):
    # Extract the "energy low [eV]" column from the dataframe and convert to float array
    nenergies = tally_dataframe["energy low [eV]"].astype(float).values
    
    # Initialize an array of energy values, starting with a small "zero" energy value
    nergs = [1e-10]
    # Append the extracted energy values to the list
    nergs.extend(nenergies.tolist())
    # Convert the list to a numpy array for numerical operations
    nergs = np.array(nergs)
    
    # Compute lethargy
    # Calculate the corrected mean values by dividing the original mean values
    # by the logarithm of the ratio of consecutive energy intervals
    mean_corrected = tally_dataframe["mean"].values / np.log(nergs[1:] / nergs[:-1])
    
    # Calculate the corrected standard deviation in a similar manner
    stddev_corrected = tally_dataframe['std. dev.'].values / np.log(nergs[1:] / nergs[:-1])
    
    # Return the corrected mean and standard deviation arrays
    return mean_corrected, stddev_corrected

def process_g_tally(tally_dataframe):
    # Extract the "energy low [eV]" column, convert to float and then to MeV
    genergies = tally_dataframe["energy low [eV]"].astype(float).values / 1e6
    
    # Initialize an array of energy values, starting with a small "zero" energy value
    gergs = [1e-10]
    # Append the extracted energy values to the list
    gergs.extend(genergies.tolist())
    # Convert the list to a numpy array for numerical operations
    gergs = np.array(gergs)
    
    # Compute per unit energy
    # Calculate the corrected mean values by dividing the original mean values
    # by the difference of consecutive energy intervals
    mean_corrected = tally_dataframe["mean"].values / (gergs[1:] - gergs[:-1])
    
    # Calculate the corrected standard deviation in a similar manner
    stddev_corrected = tally_dataframe["std. dev."].values / (gergs[1:] - gergs[:-1])
    
    # Return the corrected mean and standard deviation arrays
    return mean_corrected, stddev_corrected

def postprocess_openmc_spectra(openmc_df: pd.DataFrame) -> pd.DataFrame:
    """
    Post-process the OpenMC spectra dataframe by removing certain columns.

    Parameters:
    openmc_df (pd.DataFrame): The input dataframe from OpenMC results.

    Returns:
    pd.DataFrame: The dataframe with specified columns removed.
    """
    # Remove the specified columns if they exist
    columns_to_remove = ['surface', 'particle', 'nuclide', 'score']
    for column in columns_to_remove:
        if column in openmc_df.columns:
            del openmc_df[column]
    return openmc_df