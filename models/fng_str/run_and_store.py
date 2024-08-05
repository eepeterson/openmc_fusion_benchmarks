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
    """Run the fng_str simulation and store an hdf file in results_database/"""
    args = _parse_args()

    # run reaction rate onaxis simulation
    try:
        p = subprocess.Popen(
            ["python3", "openmc_model.py", "--reaction_rates_onaxis"])
    except:
        p = subprocess.Popen(
            ["python", "openmc_model.py", "--reaction_rates_onaxis"])
    # wait for the simulation to finish
    p.wait()

    # run reaction rate offaxis simulation
    try:
        p = subprocess.Popen(
            ["python3", "openmc_model.py", "--reaction_rates_offaxis"])
    except:
        p = subprocess.Popen(
            ["python", "openmc_model.py", "--reaction_rates_offaxis"])
    # wait for the simulation to finish
    p.wait()

    # run nuclear heating simulation
    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--heating"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--heating"])
    # wait for the simulation to finish
    p.wait()

    # read statepoint file
    onaxis_file = ofb.ResultsFromOpenmc(
        'reaction_rates_onaxis/statepoint.100.h5')
    offaxis_file = ofb.ResultsFromOpenmc(
        'reaction_rates_offaxis/statepoint.100.h5')
    heating_file = ofb.ResultsFromOpenmc('heating/statepoint.100.h5')

    # openmc hdf file
    filename = ofb.build_hdf_filename(
        'openmc', heating_file.get_openmc_version, args.xslib)

    # store activation foil results
    xaxis_name = 'Shield depth (cm)'
    xaxis_list_onaxis1 = ['0.25', '12.95', '25.95', '38.65', '43.82', '46.35',
                          '53.3', '60.05', '66.9', '73.9', '80.6', '87.25', '91.65']
    xaxis_list_onaxis2 = ['5 / 39.12', '4 / 39.12', '6 / 41.47', '3 / 41.47',
                          '8 / 41.47', '10 / 41.47', '9 / 41.47', '11 / 41.47',
                          '1 / 43.82', '7 / 43.82', '2 / 43.82']
    xaxis_list_offaxis = ['0.25', '12.95', '25.95', '38.65', '5 / 39.12', '4 / 39.12',
                          '6 / 41.47', '3 / 41.47', '8 / 41.47', '10 / 41.47', '9 / 41.47',
                          '11 / 41.47', '1 / 43.82', '7 / 43.82', '2 / 43.82']
    for foil in helpers.foil_list:
        # on axis group 1
        onaxis_file.tally_to_hdf(tally_name=f'rr_onaxis1_{foil}',
                                 normalize_over=helpers.volumes_onaxis1,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 xaxis_list=xaxis_list_onaxis1,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)

        # on axis group 2
        onaxis_file.tally_to_hdf(tally_name=f'rr_onaxis2_{foil}',
                                 normalize_over=helpers.volumes_onaxis2,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 xaxis_list=xaxis_list_onaxis2,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)

        # off axis
        offaxis_file.tally_to_hdf(tally_name=f'rr_offaxis_{foil}',
                                  normalize_over=helpers.volumes_offaxis,
                                  xs_library=args.xslib,
                                  xaxis_name=xaxis_name,
                                  xaxis_list=xaxis_list_offaxis,
                                  path_to_database='results_database',
                                  when=args.when,
                                  where=args.where)

    # store nuclear heating results
    xaxis_heating = ['46.35/SS', '53.3/SS', '60.05/SS', '66.9/SS', '73.9/SS', '80.6/SS',
                     '87.25/SS', '91.65/SS', '95.36/SS', '97.56/Cu', '99.76/SS', '101.96/Cu']
    qtld_coeff_eff3 = {'Ce': np.array([.1, .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                       'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                       'Cn': np.array([2.7, 2.78, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}
    qtld_coeff_fendl1 = {'Ce': np.array([.1, .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                         'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                         'Cn': np.array([2.7, 2.8, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}
    qtld_coeff_fendl2 = {'Ce': np.array([.1, .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                         'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                         'Cn': np.array([2.7, 2.78, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}

    # rearrange
    tally_name = 'nuclear_heating'
    openmc_df = heating_file.get_tally_dataframe(tally_name)
    tally = helpers.postprocess_openmc_heating(openmc_df, qtld_coeff_fendl2)
    d = {xaxis_name: xaxis_heating,
         'mean': tally['mean'], 'std. dev.': tally['std. dev.']}
    tally_df = pd.DataFrame(d)

    path_to_file = Path('results_database') / filename

    # write the tally in the hdf file
    tally_df.to_hdf(path_to_file, tally_name, mode='a',
                    format='table', data_columns=True, index=False)
    code_version = 'openmc-' + \
        '.'.join(map(str, heating_file.get_openmc_version))

    # write attributes to the hdf file
    with h5py.File(path_to_file, 'a') as f:
        f[tally_name + '/table'].attrs['x_axis'] = xaxis_name
        f.attrs['code_version'] = code_version
        f.attrs['xs_library'] = args.xslib.strip().replace(' ', '')
        f.attrs['batches'] = heating_file.get_batches
        f.attrs['particles_per_batch'] = heating_file.get_particles_per_batch
        f.attrs['when'] = args.when
        f.attrs['where'] = args.where


if __name__ == "__main__":
    main()
