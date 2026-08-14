import h5py
import pandas as pd
import xarray as xr
import numpy as np
import openmc_fusion_benchmarks as ofb
import os

def get_tally_dataframe(file, tally_name):

    with h5py.File(file) as f:
        df = pd.DataFrame(f[tally_name+'/table'][()]).drop(columns='index')

    return df

def save_result(df, filename, group, realization_label):
    # Convert to xarray and add dimensions
    t = xr.DataArray(
        df.values[np.newaxis, :, :],  # shape: (1, r, c)
        dims=["realization", "row", "column"],
        coords={
            "realization": [realization_label],
            "column": df.columns,
            "row": np.arange(df.shape[0]),
        },
        name=group
    )

    # Save the tally data to a netCDF file
    ofb.utils._save_result(new_result=t, filename=filename,
                    group=group, realization_label=realization_label)
    
    return 

def convert_all_files():
    folder_path = 'oktavian_al'
    # Get only file names in the directory
    file_names = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

    for file in file_names:
        if file.endswith('.h5'):
            ndf = get_tally_dataframe(os.path.join(folder_path, file), 'nspectrum')
            pdf = get_tally_dataframe(os.path.join(folder_path, file), 'gspectrum')

            new_filename = file.replace('.h5', '_csg.h5')

            save_result(ndf, new_filename, 'neutron_leakage', 'baseline')
            save_result(pdf, new_filename, 'photon_leakage', 'baseline')

# if __name__ == "__main__":
#     convert_all_files()