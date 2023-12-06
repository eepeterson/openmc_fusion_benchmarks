#!/usr/bin/env python3
import argparse
import subprocess
import openmc_sinbad_benchmarks as osb
import helpers

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

    # run simulation
    try:
        p = subprocess.Popen(["python3", "openmc_model.py"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py"])

    # wait for the simulation to finish
    p.wait()

    # read statepoint file
    openmc_file = osb.ResultsFromOpenmc('statepoint.100.h5', 'results')

    # openmc hdf file
    filename = 'openmc_' + args.xslib + '.h5'

    # store activation foil results
    x_axis = 'Detector No.'
    for foil in helpers.foil_list:
        openmc_file.tally_to_hdf(hdf_file=filename, tally_name=f'rr_{foil}',
                                    normalize_over=helpers.detector_volume,
                                    xs_library=args.xlib, x_axis=x_axis,
                                    path_to_database='results_database',
                                    when=args.when,
                                    where=args.where)
        
    # store spectrometer results
    x_axis = 'Detector No.'
    for dp in helpers.detector_list:
        openmc_file.tally_to_hdf(hdf_file=filename, tally_name=f'nspectrum_detector{dp}',
                                    normalize_over=helpers.detector_volume,
                                    xs_library=args.xlib, x_axis=x_axis,
                                    path_to_database='results_database',
                                    when=args.when,
                                    where=args.where)

if __name__ == "__main__":
    main()