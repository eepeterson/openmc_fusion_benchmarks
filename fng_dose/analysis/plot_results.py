from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
from uncertainties import ufloat

import openmc.deplete

mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['font.size'] = 12.0

# Experimental data
exp_data = np.loadtxt('experiment.txt')
exp_data[:, 0] *= 365.25

# Load OpenMC results
results_dir = Path('../results/results_campaign1_endf80_icrp74')
n = 18  # number of cooling times
results = openmc.deplete.Results(results_dir / 'depletion_results.h5')
times = results.get_times()
times = times[-n:] - times[-n-1]
dose_openmc = np.load(results_dir / 'dose.npy')

# Plot against experiment
fig, ax = plt.subplots()
ax.loglog(exp_data[:, 0], exp_data[:, 1], label='Experiment')
ax.loglog(times, dose_openmc[:, 0], 'kx', label='OpenMC')
ax.set_xlabel('Time [d]')
ax.set_ylabel('Dose rate [Sv/h]')
ax.grid(which='both')
ax.legend(framealpha=1.0)
plt.savefig('dose_campaign1.pdf', bbox_inches='tight')
#plt.show()

# Compare at selected 5 points vs experiment and results from Eade
dose_experiment = np.array([2.46e-6, 6.99e-7, 4.95e-7, 4.16e-7, 3.16e-7])
mcr2s_2010 = [0.883, 1.115, 1.230, 1.233, 1.214]
mcr2s_2020 = [2.95e-6, 8.92e-7, 7.18e-7, 6.18e-7, 4.68e-7]
comparison_index = [5, 10, 13, 16, 17]
selected = dose_openmc[comparison_index, 0]
for mean, unc in dose_openmc[comparison_index]:
    print(ufloat(mean, unc))

fig, ax = plt.subplots()
plt.plot(range(5), mcr2s_2010, 'rs', label='MCR2S (2010)')
plt.plot(range(5), mcr2s_2020/dose_experiment, 'ro', label='MCR2S (2020)')
plt.plot(range(5), selected/dose_experiment, 'bx', label='OpenMC')
xmin, xmax = -0.3, 4.3
plt.plot([xmin, xmax], [1., 1.], 'k--')
sigma = 0.1
verts = [(xmin, 1. - sigma), (xmin, 1. + sigma), (xmax, 1. + sigma), (xmax, 1. - sigma)]
poly = Polygon(verts, facecolor='gray', alpha=0.5)
ax.add_patch(poly)
ax.set_xlabel('Time', fontsize=14)
ax.set_ylabel('C/E', fontsize=14)
#ax.set_ylim(0.7, 1.3)
ax.set_xticks(range(5), ['1 d', '7 d', '15 d', '30 d', '60 d'])
ax.set_xlim(xmin, xmax)
ax.legend(loc='upper left', framealpha=1.0)
ax.grid(which='both')
plt.savefig('dose_campaign1_ce.pdf', bbox_inches='tight')
#plt.show()
