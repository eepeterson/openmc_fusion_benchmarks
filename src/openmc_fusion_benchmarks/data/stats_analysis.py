import re
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy.stats import norm

import openmc_fusion_benchmarks as ofb


class ResultsFromTMC(ofb.ResultsFromDatabase):
    def __init__(self, file):
        super().__init__(file)

        # self.filename = file.strip().split('/')
        # self.filepath = Path(file)

    # def list_tallies(self):
    #     # Create a set to store unique base names
    #     base_names = set()

    #     # Regular expression to capture the pattern 'name_N'
    #     pattern = re.compile(r"(.+?)_(\d+)$")

    #     # Open the HDF5 file
    #     with h5py.File(self.filepath, 'r') as hdf:
    #         # List all datasets in the HDF5 file
    #         dataset_names = list(hdf.keys())

    #         # Loop through the dataset names
    #         for name in dataset_names:
    #             match = pattern.match(name)
    #             if match:
    #                 # Capture the base name
    #                 base_name = match.group(1)
    #                 # Add to the set of unique base names
    #                 base_names.add(base_name)

    def find_nsamples(self, tally):

        # Open the HDF5 file
        with h5py.File(self.filepath, 'r') as hdf:
            # List all datasets in the HDF5 file
            dataset_names = list(hdf.keys())
            # Initialize to track the highest N
            max_n = -1
            # Loop through the datasets to find the ones that match the pattern base_name_N
            for name in dataset_names:
                if name.startswith(tally + '_'):
                    # Extract N and update max_n if N is larger
                    n = int(name.split('_')[-1])
                    if n > max_n:
                        max_n = n

        if max_n == -1:
            msg = f"No datasets found with base name '{tally}' in the HDF5 file."
            raise ValueError(msg)
        else:
            return max_n

    def get_means(self, tally):
        mean = []
        nsamples = self.find_nsamples(tally)
        with h5py.File(self.filepath, 'r') as f:
            for n in range(nsamples):
                df = pd.DataFrame(
                    f[tally+f'_{n}'+'/table'][()]).drop(columns='index')
                mean.append(df['mean'][0])

        return mean

    def get_std_devs(self, tally):
        std_dev = []
        nsamples = self.find_nsamples(tally)
        with h5py.File(self.filepath, 'r') as f:
            for n in range(nsamples):
                df = pd.DataFrame(
                    f[tally+f'_{n}'+'/table'][()]).drop(columns='index')
                std_dev.append(df['std. dev.'][0])

        return std_dev


def get_gauss(mu, sigma):
    x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
    y = norm.pdf(x, mu, sigma)
    y = y / max(y)

    return x, y
