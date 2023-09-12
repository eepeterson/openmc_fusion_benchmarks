from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import openmc.deplete
from uncertainties import ufloat

mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 12.0
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['savefig.bbox'] = 'tight'

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
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 7))
ax1.loglog(exp_data[:, 0], exp_data[:, 1], 'k-', label='Experiment', lw=0.75)
ax1.loglog(times, dose_openmc[:, 0], 'bx', label='OpenMC')
ax1.set_ylabel('Dose rate [Sv/h]')
ax1.grid(which='both')
ax1.legend(framealpha=1.0)

# Compare at selected 5 points vs experiment and results from Eade
dose_experiment = np.array([
    1.30e-4, # 1 hour = 1.14e-4 a
    3.55e-5, # 6 hours = 6.85e-4 a
    9.19e-6, # 12 hours = 1.37e-3 a
    4.77e-6, # 16 hours = 1.83e-3 a
    3.23e-6, # 20 hours = 2.28e-3 a
    2.46e-6, # 1 day = 2.74e-3 a
    1.46e-6, # 2 days = 5.48e-3 a
    1.18e-6, # 3 days = 8.22-3 a
    9.68e-7, # 4 days = 1.10e-2 a
    8.26e-7, # 5 days = 1.37e-2 a
    6.99e-7, # 7 days = 1.92e-2 a
    5.96e-7, # 9 days = 2.47e-2 a
    5.09e-7, # 12 days = 3.29e-2 a
    4.95e-7, # 15 days = 4.11e-2 a
    4.85e-7, # 18 days = 4.93e-2 a
    4.55e-7, # 21 days = 5.75e-2 a
    4.16e-7, # 30 days = 8.22e-2 a
    3.16e-7, # 60 days = 1.64e-1 a
])

dose_selected = np.array([2.46e-6, 6.99e-7, 4.95e-7, 4.16e-7, 3.16e-7])
comparison_index = [5, 10, 13, 16, 17]
mcr2s_times = [times[i] for i in comparison_index]
mcr2s_2010 = [0.883, 1.115, 1.230, 1.233, 1.214]
mcr2s_2020 = [2.95e-6, 8.92e-7, 7.18e-7, 6.18e-7, 4.68e-7]
selected = dose_openmc[comparison_index, 0]
for mean, unc in dose_openmc[comparison_index]:
    print(ufloat(mean, unc))

ax2.plot(mcr2s_times, mcr2s_2010, 'ro', fillstyle='none', label='MCR2S (2010)')
ax2.plot(mcr2s_times, mcr2s_2020/dose_selected, 'r^', label='MCR2S (2020)')
ax2.semilogx(times, dose_openmc[:, 0]/dose_experiment, 'bx', label='OpenMC')
xmin, xmax = exp_data[:, 0].min(), exp_data[:, 0].max()
ax2.plot([xmin, xmax], [1., 1.], 'k--', lw=0.75)
sigma = 0.1
verts = [(xmin, 1. - sigma), (xmin, 1. + sigma), (xmax, 1. + sigma), (xmax, 1. - sigma)]
poly = Polygon(verts, facecolor='gray', alpha=0.5)
ax2.add_patch(poly)
ax2.set_xlabel('Time [d]', fontsize=14)
ax2.set_ylabel('C/E', fontsize=14)
ax2.set_ylim(0.8, 1.5)
ax2.legend(loc='upper left', framealpha=1.0)
ax2.grid(which='both')
fig.savefig('dose_campaign1.pdf')
