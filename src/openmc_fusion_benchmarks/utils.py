import numpy as np
import pandas as pd
from typing import Iterable


def rescale_to_lethargy(df: pd.DataFrame) -> pd.DataFrame:
    """Rescale the mean and std. dev. of a DataFrame to lethargy values.
    Useful when it comes to plot a particle energy spectrum in lethargy scale.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing an energy spectrum tally. it needs to contain
        'mean' and 'std. dev.' columns, as well as 'energy low [eV]' and 
        'energy high [eV]' columns.

    Returns
    -------
    pd.DataFrame
        DataFrame with the mean and std. dev. rescaled to lethargy values.
    """
    # Copy the DataFrame to avoid modifying the original
    df_lethargy = df.copy()

    # Calculate the lethargy of the energy bins
    lethargy = np.log(
        df_lethargy['energy high [eV]'] / df_lethargy['energy low [eV]'])

    # Rescale the mean and std. dev. to lethargy
    df_lethargy['mean'] = df_lethargy['mean'] / lethargy
    df_lethargy['std. dev.'] = df_lethargy['std. dev.'] / lethargy

    return df_lethargy


def rebin_spectrum(df: pd.DataFrame, energy_low: Iterable, energy_high: Iterable) -> pd.DataFrame:
    """Rebin an energy spectrum tally to new energy bins. The new energy bins
    are defined by the energy_low and energy_high arrays. The function
    calculates the mean and std. dev. of the energy spectrum in the new bins.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing an energy spectrum tally. It needs to contain
        'mean' and 'std. dev.' columns, as well as 'energy low [eV]' and 
        'energy high [eV]' columns.
    energy_low : Iterable
        Iterable containing the lower bounds of the new energy bins.
    energy_high : Iterable
        Iterable containing the upper bounds of the new energy bins.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the rebinned energy spectrum tally. It has the
        same columns as the input DataFrame, but the mean and std. dev. are
        calculated in the new energy bins.
    """

    # Ensure that energy_low and energy_high are numpy arrays
    energy_low = np.array(energy_low)
    energy_high = np.array(energy_high)

    # Create the target DataFrame from the energy bounds
    df_to = pd.DataFrame(
        {'energy low [eV]': energy_low, 'energy high [eV]': energy_high})

    # Prepare arrays to store the rebinned mean and std. dev.
    rebinned_probs = np.zeros(len(df_to))
    rebinned_std_devs = np.zeros(len(df_to))

    for i in range(len(df_to)):
        bin_low = df_to['energy low [eV]'].iloc[i]
        bin_high = df_to['energy high [eV]'].iloc[i]

        # Find overlapping bins
        mask = (df['energy high [eV]'] > bin_low) & (
            df['energy low [eV]'] < bin_high)

        if mask.any():
            # Calculate overlapping range for each bin
            overlaps = df[mask].copy()
            overlaps['overlap_low'] = np.maximum(
                overlaps['energy low [eV]'], bin_low)
            overlaps['overlap_high'] = np.minimum(
                overlaps['energy high [eV]'], bin_high)
            overlaps['overlap_width'] = overlaps['overlap_high'] - \
                overlaps['overlap_low']

            # Calculate the weighted mean and std. dev. for the new bin
            weighted_probs = overlaps['mean'] * overlaps['overlap_width']
            rebinned_probs[i] = weighted_probs.sum() / (bin_high - bin_low)

            # Calculate the weighted std. dev. for the new bin
            weighted_vars = (overlaps['std. dev.']**2) * \
                (overlaps['overlap_width']**2)
            rebinned_std_devs[i] = np.sqrt(
                weighted_vars.sum()) / (bin_high - bin_low)

    # Create a new DataFrame with the rebinned mean and std. dev.
    df_rebinned = df_to.copy()
    df_rebinned['mean'] = rebinned_probs
    df_rebinned['std. dev.'] = rebinned_std_devs
    return df_rebinned


def get_nonzero_energy_interval(df: pd.DataFrame) -> tuple:
    """Function that takes in a pandas dataframe with a energy spectrum tally
    and returns the energy interval where the mean is nonzero. Useful to set
    the correct bounds for plotting the energy spectrum.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe containing an energy spectrum tally.

    Returns
    -------
    tuple
        tuple with the energy low and high values of the interval where the
        mean is nonzero. If the mean is zero for all energy bins, the function
        returns the energy low and high values of the first and last energy
        bins.
    """
    # Identify the first and last index where 'mean' is nonzero
    nonzero_indices = df.index[df['mean'] != 0].tolist()

    if not nonzero_indices:
        # If there are no nonzero values in 'mean'
        first_energy_low = np.array(df['energy low [eV]'])[0]
        last_energy_high = np.array(df['energy high [eV]'])[-1]
        return first_energy_low, last_energy_high

    first_nonzero_idx = nonzero_indices[0]
    last_nonzero_idx = nonzero_indices[-1]

    # Get the corresponding energy low and high values
    first_energy_low = df.at[first_nonzero_idx, 'energy low [eV]']
    last_energy_high = df.at[last_nonzero_idx, 'energy high [eV]']

    return first_energy_low, last_energy_high
