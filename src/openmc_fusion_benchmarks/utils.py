import numpy as np
import pandas as pd
from typing import Iterable


def rebin_spectrum(df: pd.DataFrame, energy_low: Iterable, energy_high: Iterable):
    # Ensure that energy_low and energy_high are numpy arrays
    energy_low = np.array(energy_low)
    energy_high = np.array(energy_high)

    # Create the target DataFrame from the energy bounds
    df_to = pd.DataFrame(
        {'energy low [eV]': energy_low, 'energy high [eV]': energy_high})

    # Prepare arrays to store the rebinned probabilities and their standard deviations
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

            # Calculate the weighted probability and standard deviation for the new bin
            weighted_probs = overlaps['mean'] * overlaps['overlap_width']
            rebinned_probs[i] = weighted_probs.sum() / (bin_high - bin_low)

            # Calculate the weighted standard deviation for the new bin
            weighted_vars = (overlaps['std. dev.']**2) * \
                (overlaps['overlap_width']**2)
            rebinned_std_devs[i] = np.sqrt(
                weighted_vars.sum()) / (bin_high - bin_low)

    # Create a new DataFrame with the rebinned probabilities and standard deviations
    df_rebinned = df_to.copy()
    df_rebinned['mean'] = rebinned_probs
    df_rebinned['std. dev.'] = rebinned_std_devs
    return df_rebinned
