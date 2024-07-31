#!/usr/bin/env python3
import argparse
import subprocess
import openmc_fusion_benchmarks as ofb
from pathlib import Path
import helpers
import h5py

# ignore NaturalNameWarnings
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xslib", type=str)
    parser.add_argument("-t", "--when", type=str, default='n/a')
    parser.add_argument("-w", "--where", type=str, default='n/a')

    args = parser.parse_args()

    return args


def main():
    """Run the fns_duct simulation and store an hdf file in results_database/"""
    args = _parse_args()

    if args.xslib is None:
        msg = """Please enter the used nuclear data library name to be used in the simulation.
            Specify it with the argument -x or --xslib"""
        raise ValueError(msg)

    # run simulation
    try:
        p = subprocess.Popen(["python3", "openmc_model.py"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py"])

    # wait for the simulation to finish
    p.wait()

    # read statepoint file
    openmc_file = ofb.ResultsFromOpenmc('statepoint.100.h5', 'results')

    # generate openmc hdf file
    filename = ofb.build_hdf_filename(
        'openmc', openmc_file.get_openmc_version, args.xslib)
    
    # store nuclear spectra results
    # # rearrange
    xaxis_name = 'energy low [eV]'
    n_tally_name = 'nspectrum'
    n_openmc_df = openmc_file.get_tally_dataframe(n_tally_name)
    n_tally_df = helpers.postprocess_openmc_spectra(n_openmc_df)

    path_to_file = Path('results_database') / filename

    # write the tally in the hdf file
    n_tally_df.to_hdf(path_to_file, n_tally_name, mode='a',
                    format='table', data_columns=True, index=False)
    code_version = 'openmc-' + \
        '.'.join(map(str, openmc_file.get_openmc_version))

    # write attributes to the hdf file
    with h5py.File(path_to_file, 'a') as f:
        f[n_tally_name + '/table'].attrs['x_axis'] = xaxis_name
        f.attrs['code_version'] = code_version
        f.attrs['xs_library'] = args.xslib.strip().replace(' ', '')
        f.attrs['batches'] = openmc_file.get_batches
        f.attrs['particles_per_batch'] = openmc_file.get_particles_per_batch
        f.attrs['when'] = args.when
        f.attrs['where'] = args.where

    xaxis_name = 'energy low [eV]'
    g_tally_name = 'gspectrum'
    g_openmc_df = openmc_file.get_tally_dataframe(g_tally_name)
    g_tally_df = helpers.postprocess_openmc_spectra(g_openmc_df)

    path_to_file = Path('results_database') / filename

    # write the tally in the hdf file
    g_tally_df.to_hdf(path_to_file, g_tally_name, mode='a',
                    format='table', data_columns=True, index=False)
    code_version = 'openmc-' + \
        '.'.join(map(str, openmc_file.get_openmc_version))

    # write attributes to the hdf file
    with h5py.File(path_to_file, 'a') as f:
        f[g_tally_name + '/table'].attrs['x_axis'] = xaxis_name
        f.attrs['code_version'] = code_version
        f.attrs['xs_library'] = args.xslib.strip().replace(' ', '')
        f.attrs['batches'] = openmc_file.get_batches
        f.attrs['particles_per_batch'] = openmc_file.get_particles_per_batch
        f.attrs['when'] = args.when
        f.attrs['where'] = args.where
        
if __name__ == "__main__":
    main()
