import argparse
from math import pi

import numpy as np
import openmc
from prettytable import PrettyTable
from uncertainties import ufloat


parser = argparse.ArgumentParser()
parser.add_argument('statepoint')
args = parser.parse_args()

with openmc.StatePoint(args.statepoint) as sp:
    ni_rr_tally = list(sp.tallies.values())[0]

# Calculate volume for each foil
r = 0.9
thickness = np.array([0.2, 0.2, 0.1, 0.2, 0.2, 0.2])
volume = pi*r**2 * thickness

n2n_mean = ni_rr_tally.get_values(scores=['(n,2n)']).flatten() / volume
np_mean = ni_rr_tally.get_values(scores=['(n,p)']).flatten() / volume
n2n_std = ni_rr_tally.get_values(scores=['(n,2n)'], value='std_dev').flatten() / volume
np_std = ni_rr_tally.get_values(scores=['(n,p)'], value='std_dev').flatten() / volume

# Experimental values for reaction rates
experiment_np_rate = np.array([21.9, 5.19, 4.13, 8.48, 7.86, 5.15])*1e-6
experiment_n2n_rate = np.array([28.4, 3.94, 2.07, 4.92, 4.71, 3.64])*1e-7

table = PrettyTable(["Foil", "Experiment", "OpenMC"])
for i, (exp, mean, std) in enumerate(zip(experiment_np_rate, np_mean, np_std)):
    table.add_row([i+1, f'{exp:.2e}', ufloat(mean, std)])
print("(n,p) reaction rate")
print(table)

table = PrettyTable(["Foil", "Experiment", "OpenMC"])
for i, (exp, mean, std) in enumerate(zip(experiment_n2n_rate, n2n_mean, n2n_std)):
    table.add_row([i+1, f'{exp:.2e}', ufloat(mean, std)])
print("\n(n,2n) reaction rate")
print(table)
