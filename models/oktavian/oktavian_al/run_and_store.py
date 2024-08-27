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
    parser.add_argument("-x", "--xslib", type=str,
                        help="String with Cross section library name and version (e.g. 'FENDL-2.3')")
    parser.add_argument("-t", "--when", type=str, default='n/a',
                        help="String with the month and year the simulation is run as (e.g. 'June 2021')")
    parser.add_argument("-w", "--where", type=str, default='n/a',
                        help="String with the place/institution where the simulation is run (e.g. 'MIT-PSFC')")

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

    xaxis_name = 'energy low [eV]'
    openmc_file.tally_to_hdf(tally_name=f'nspectrum',
                             normalize_over=helpers.surface,
                             xs_library=args.xslib,
                             xaxis_name=xaxis_name,
                             path_to_database='results_database',
                             when=args.when,
                             where=args.where)
    
    openmc_file.tally_to_hdf(tally_name=f'gspectrum',
                             normalize_over=helpers.surface,
                             xs_library=args.xslib,
                             xaxis_name=xaxis_name,
                             path_to_database='results_database',
                             when=args.when,
                             where=args.where)
        
if __name__ == "__main__":
    main()
