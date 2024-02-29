#!/usr/bin/env python3
import argparse
import subprocess
import openmc_fusion_benchmarks as ofb
import helpers

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

    # store activation foil results
    xaxis_name = 'Detector No.'
    xaxis_list = ['1',  '2',  '3',  '4',  '5',
                  '6',  '7',  '8',  '9', '10', '11']
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
    for dp in helpers.detector_list:
        openmc_file.tally_to_hdf(tally_name=f'nspectrum_{dp}',
                                 normalize_over=helpers.detector_volume,
                                 xs_library=args.xslib,
                                 xaxis_name=xaxis_name,
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)


if __name__ == "__main__":
    main()
