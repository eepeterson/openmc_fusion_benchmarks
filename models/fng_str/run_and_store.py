#!/usr/bin/env python3
import argparse
import subprocess
import openmc_sinbad_benchmarks as osb
import helpers
import numpy as np
import pandas as pd
from pathlib import Path
import h5py

# ignore NaturalNameWarnings
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)

_xaxis_heating = ['46.35/SS', '53.3/SS', '60.05/SS', '66.9/SS', '73.9/SS', '80.6/SS', 
          '87.25/SS', '91.65/SS', '95.36/SS', '97.56/Cu', '99.76/SS', '101.96/Cu']
_qtld_coeff_eff3 = {'Ce': np.array([.1 , .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                    'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                    'Cn': np.array([2.7, 2.78, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}
_qtld_coeff_fendl1 = {'Ce': np.array([.1 , .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                    'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                    'Cn': np.array([2.7, 2.8, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}
_qtld_coeff_fendl2 = {'Ce': np.array([.1 , .07, .05, .05, .05, .05, .05, .05, .05, .05, .05, .05]),
                    'Cp': np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05]),
                    'Cn': np.array([2.7, 2.78, 2.9, 3., 3.09, 3.19, 3.29, 3.39, 3.5, 4.54, 3.65, 5.12])}

def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xslib", type=str)
    parser.add_argument("-t", "--when", type=str, default='n/a')
    parser.add_argument("-w", "--where", type=str, default='n/a')

    args = parser.parse_args()

    return args

def main():
    """Run the fng_str simulation and store an hdf file in results_database/"""
    args = _parse_args()

    # run reaction rate onaxis simulation
    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--reaction_rates_onaxis"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--reaction_rates_onaxis"])
    # wait for the simulation to finish
    p.wait()

    # run reaction rate offaxis simulation
    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--reaction_rates_offaxis"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--reaction_rates_offaxis"])
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
    onaxis_file = osb.ResultsFromOpenmc('statepoint.100.h5', 'reaction_rates_onaxis')
    offaxis_file = osb.ResultsFromOpenmc('statepoint.100.h5', 'reaction_rates_offaxis')
    heating_file = osb.ResultsFromOpenmc('statepoint.100.h5', 'heating')

    # openmc hdf file
    filename = osb.build_hdf_filename('openmc', args.xslib)

    # store activation foil results
    x_axis = 'Shield depth (cm)'
    for foil in helpers.foil_list:
        # on axis group 1
        onaxis_file.tally_to_hdf(hdf_file=filename, tally_name=f'rr_onaxis1_{foil}',
                                    normalize_over=helpers.volumes_onaxis1,
                                    xs_library=args.xlib, x_axis=x_axis,
                                    path_to_database='results_database',
                                    when=args.when,
                                    where=args.where)
        
        # on axis group 2
        onaxis_file.tally_to_hdf(hdf_file=filename, tally_name=f'rr_onaxis2_{foil}',
                                    normalize_over=helpers.volumes_onaxis2,
                                    xs_library=args.xlib, x_axis=x_axis,
                                    path_to_database='results_database',
                                    when=args.when,
                                    where=args.where)
        
        # off axis
        offaxis_file.tally_to_hdf(hdf_file=filename, tally_name=f'rr_offaxis_{foil}',
                                    normalize_over=helpers.volumes_offaxis,
                                    xs_library=args.xlib, x_axis=x_axis,
                                    path_to_database='results_database',
                                    when=args.when,
                                    where=args.where)
        
    # store nuclear heating results
    # rearrange
    tally_name = 'nuclear_heating'
    openmc_df = heating_file.get_tally_dataframe(tally_name)
    tally = helpers.postprocess_openmc_heating(openmc_df, _qtld_coeff_fendl2)
    d = {x_axis:_xaxis_heating, 'mean':tally['mean'], 'std. dev.':tally['std. dev.']}
    tally_df = pd.Dataframe(d)

    path_to_file = Path('results_database') / filename

    # write the tally in the hdf file
    tally_df.to_hdf(path_to_file, tally_name, mode='a',
                    format='table', data_columns=True, index=False)
    code_version = 'openmc-' + '.'.join(map(str, heating_file.get_openmc_version))

    # write attributes to the hdf file
    with h5py.File(path_to_file, 'a') as f:
        f[tally_name + '/table'].attrs['x_axis'] = x_axis
        f.attrs['code_version'] = code_version
        f.attrs['xs_library'] = args.xlib.strip.replace(' ','')
        f.attrs['batches'] = heating_file.get_batches
        f.attrs['particles_per_batch'] = heating_file.get_particles_per_batch
        f.attrs['when'] = args.when
        f.attrs['where'] = args.where


if __name__ == "__main__":
    main()