import argparse
from math import pi

from matplotlib.patches import Polygon
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import openmc
from prettytable import PrettyTable
from uncertainties import ufloat

mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['font.size'] = 12.0

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

ce_n2n = n2n_mean/experiment_n2n_rate
ce_n2n_unc = n2n_std/experiment_n2n_rate
ce_np = np_mean/experiment_np_rate
ce_np_unc = np_std/experiment_np_rate

fig, ax = plt.subplots()
ax.errorbar([1, 2, 3, 4, 5, 6], ce_np, ce_np_unc, ls='none', marker='x', label='(n,p)')
ax.errorbar([1, 2, 3, 4, 5, 6], ce_n2n, ce_n2n_unc, ls='none', marker='o', label='(n,2n)')
poly = Polygon([(0, 0.95), (0, 1.05), (7, 1.05), (7, 0.95)], facecolor='gray', edgecolor=None, alpha=0.2)
ax.add_patch(poly)
ax.legend(loc='upper right', framealpha=1.0)
ax.set_xlim(0.5, 6.5)
ax.set_ylim(0.75, 1.15)
ax.set_xlabel('Foil')
ax.set_ylabel('C/E')
ax.grid(True, axis='y')
plt.savefig('ni_foil_ce.pdf')
