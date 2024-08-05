#!/usr/bin/env python3
import argparse
import subprocess
import openmc_fusion_benchmarks as ofb
import helpers
import numpy as np
import pandas as pd
from pathlib import Path
import h5py

# ignore NaturalNameWarnings
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xslib", type=str,
                        help="Strign with Cross section library name and version (e.g. 'FENDL-2.3')")
    parser.add_argument("-t", "--when", type=str, default='n/a',
                        help="String with the month and year the simulation was run as (e.g. 'June 2021')")
    parser.add_argument("-w", "--where", type=str, default='n/a',
                        help="String with the place/institution where the simulation is run (e.g. 'MIT-PSFC')")

    args = parser.parse_args()

    return args


def main():
    """Run the fns_clean_w simulation and store an hdf file in results_database/"""
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
    openmc_file = ofb.ResultsFromOpenmc('results/statepoint.100.h5')

    # openmc hdf file
    filename = ofb.build_hdf_filename(
        'openmc', openmc_file.get_openmc_version, args.xslib)

    # store activation foil results
    xaxis_name = 'Shield depth (cm)'
    xaxis_list = ['0.0', '7.6', '22.8', '38.0', '50.7']
    for foil in helpers.foil_list:
        openmc_file.tally_to_hdf(tally_name=f'rr_{foil}',
                                 normalize_over=helpers.detector_volume,
                                 xs_library=args.xslib, xaxis_name=xaxis_name,
                                 xaxis_list=xaxis_list,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)

    # store spectrometer results
    xaxis_name = 'Energy low [eV]'
    for dp, v in zip(helpers.detector_list, helpers.detector_volume[1:-1]):
        # ne213 neutron spectrometer
        openmc_file.tally_to_hdf(tally_name=f'nspectrum_ne213_{dp}',
                                 normalize_over=v,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)
        # prc neutron spectrometer
        openmc_file.tally_to_hdf(tally_name=f'nspectrum_prc_{dp}',
                                 normalize_over=v,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)
        # bc537 gamma spectrometer
        openmc_file.tally_to_hdf(tally_name=f'gspectrum_bc537_{dp}',
                                 normalize_over=v,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)

    # rearrange
    tally_name = 'nuclear_heating'
    openmc_df = openmc_file.get_tally_dataframe(tally_name)
    tally = helpers.postprocess_openmc_heating(openmc_df)
    d = {xaxis_name: xaxis_list[:-1],
         'mean': tally['mean'], 'std. dev.': tally['std. dev.']}
    tally_df = pd.DataFrame(d)

    path_to_file = Path('results_database') / filename

    # write the tally in the hdf file
    tally_df.to_hdf(path_to_file, tally_name, mode='a',
                    format='table', data_columns=True, index=False)
    code_version = 'openmc-' + \
        '.'.join(map(str, openmc_file.get_openmc_version))

    # write attributes to the hdf file
    with h5py.File(path_to_file, 'a') as f:
        f[tally_name + '/table'].attrs['x_axis'] = xaxis_name
        f.attrs['code_version'] = code_version
        f.attrs['xs_library'] = args.xslib.strip().replace(' ', '')
        f.attrs['batches'] = openmc_file.get_batches
        f.attrs['particles_per_batch'] = openmc_file.get_particles_per_batch
        f.attrs['when'] = args.when
        f.attrs['where'] = args.where


if __name__ == "__main__":
    main()
