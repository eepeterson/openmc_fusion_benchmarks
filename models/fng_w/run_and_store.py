#!/usr/bin/env python3
import argparse
import subprocess
import openmc_fusion_benchmarks as ofb
import helpers
import pandas as pd
from pathlib import Path

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
                        help="String with the month and year the simulation is run as (e.g. 'June 2021')")
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
            ["python3", "openmc_model.py", "--reaction_rates"])
    except:
        p = subprocess.Popen(
            ["python", "openmc_model.py", "--reaction_rates"])
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
    reaction_rates_file = ofb.ResultsFromOpenmc(
        'reaction_rates/statepoint.100.h5')
    heating_file = ofb.ResultsFromOpenmc('heating/statepoint.100.h5')

    # generate openmc hdf file
    filename = ofb.build_hdf_filename(
        'openmc', heating_file.get_openmc_version, args.xslib)

    # store activation foil results
    xaxis_name = 'Shield depth (cm)'
    for foil in helpers.foil_list:
        reaction_rates_file.tally_to_hdf(tally_name=f'rr_{foil}',
                                         normalize_over=helpers.foil_volumes[foil],
                                         xs_library=args.xslib,
                                         xaxis_name=xaxis_name,
                                         xaxis_list=helpers.xaxis_rr[foil],
                                         path_to_database='results_database',
                                         when=args.when,
                                         where=args.where)

    # store nuclear heating results
    # # rearrange
    tally_name = 'nuclear_heating'
    openmc_df = heating_file.get_tally_dataframe(tally_name)
    tally = helpers.postprocess_openmc_heating(openmc_df)
    d = {xaxis_name: helpers.xaxis_heating,
         'mean': tally['mean'], 'std. dev.': tally['std. dev.']}
    tally_df = pd.DataFrame(d)

    path_to_file = Path('results_database') / filename
    code_version = 'openmc-' + \
        '.'.join(map(str, heating_file.get_openmc_version))

    xs_library = args.xslib.strip().replace(' ', '')
    ofb.to_hdf(tally_df, path_to_file, tally_name, xs_library, xaxis_name,
               args.when, args.where, code_version,
               heating_file.get_batches, heating_file.get_particles_per_batch)


if __name__ == "__main__":
    main()
