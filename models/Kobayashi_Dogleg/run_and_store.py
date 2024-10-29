#!/usr/bin/env python3
import argparse
import subprocess
import openmc_fusion_benchmarks as ofb
import helpers
import numpy as np
import h5py
from pathlib import Path


# ignore NaturalNameWarnings
import warnings
from tables import NaturalNameWarning
warnings.filterwarnings('ignore', category=NaturalNameWarning)


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--when", type=str, default='n/a')
    parser.add_argument("-w", "--where", type=str, default='n/a')
    parser.add_argument("-x", "--xslib", type=str, default='mgxs')

    args = parser.parse_args()

    return args


def main():
    """Run the Kobayashi Dogleg simulation and store an hdf file in results_database/"""
    args = _parse_args()
    
    # generate weight windows tallies
    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--Problem_I","--weight_windows"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--Problem_I","--weight_windows"])

    # wait for the simulation to finish
    p.wait()

    # run simulation without scattering
    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--Problem_I"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--Problem_I"])

    # wait for the simulation to finish
    p.wait()

    try:
        p = subprocess.Popen(["python3", "openmc_model.py", "--weight_windows"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py", "--weight_windows"])

    # wait for the simulation to finish
    p.wait()

    # run simulation with scattering
    try:
        p = subprocess.Popen(["python3", "openmc_model.py"])
    except:
        p = subprocess.Popen(["python", "openmc_model.py"])

    # wait for the simulation to finish
    p.wait()
    
    print("Reading Statepoint file")
    # read statepoint file
    ProblemI_file = ofb.ResultsFromOpenmc('statepoint.20.h5', 'Problem_I')
    # read statepoint file
    ProblemII_file = ofb.ResultsFromOpenmc('statepoint.20.h5', 'Problem_II')

    # openmc hdf file
    filename = ofb.build_hdf_filename(
        'openmc', ProblemII_file.get_openmc_version, 'mgxs')

    # store case results
    xaxis_name = "Sample Point"
    #X-axis for case 3A
    xaxis_list = [[str(x) for x in np.linspace(5,95,10)]]
    #X-axis for case 3B
    xaxis_list.append([str(x) for x in np.linspace(5,65,6)])
    #X-axis for case 3C
    xaxis_list.append([str(x) for x in np.linspace(5,65,6)])
    
    for i,case in enumerate(helpers.case_list):
        print("Output to hdf ",case)
        ProblemI_file.tally_to_hdf(tally_name="Problem_I_"+case,
                                 normalize_over=helpers.detector_volume,
                                 xs_library="mgxs", xaxis_name=xaxis_name,
                                 xaxis_list=xaxis_list[i],
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)

        ProblemII_file.tally_to_hdf(tally_name="Problem_II_"+case,
                                 normalize_over=helpers.detector_volume,
                                 xs_library="mgxs", xaxis_name=xaxis_name,
                                 xaxis_list=xaxis_list[i],
                                 path_to_database='results_database',
                                 when=args.when,
                                 where=args.where)


        
        path_to_file = Path('results_database') / filename
        code_version = 'openmc-' + \
            '.'.join(map(str, ProblemII_file.get_openmc_version))
        
        # write attributes to the hdf file
        with h5py.File(path_to_file, 'a') as f:
            f[case + '/table'].attrs['x_axis'] = xaxis_name
            f.attrs['code_version'] = code_version
            f.attrs['xs_library'] = args.xslib.strip().replace(' ', '')
            f.attrs['batches'] = ProblemII_file.get_batches
            f.attrs['particles_per_batch'] = ProblemII_file.get_particles_per_batch
            f.attrs['when'] = args.when
            f.attrs['where'] = args.where
        
        
if __name__ == "__main__":
    main()
