import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Experimental data
exp_data = np.loadtxt('experiment.txt')
exp_data[:, 0] *= 365.25

# Calculated data
times = [1, 7, 15, 30, 60]
#sv = [1.0e-6, 3.91e-7, 2.91e-7, 2.87e-7, 1.87e-7]
#sv = [1.862e-6, 1.482e-6, 4.53e-7, 3.07e-7, 4.02e-7]
#sv = [1.843e-6, 1.519e-6, 4.35e-7, 3.001e-7, 4.029e-7]
#sv = [9.53e-7, 3.544e-7, 3.0e-7, 2.127e-7, 1.581e-7]
tendl = [1.201e-6, 4.30e-7, 3.633e-7, 2.615e-7, 1.933e-7]
endf8 = [1.599e-6, 8.245e-7, 4.568e-7, 3.044e-7, 1.803e-7]
#endf8_fispact = [9.370e-7, 3.188e-7, 2.443e-7, 2.060e-7, 1.598e-7]
endf7 = [9.86e-7, 4.35e-7, 2.83e-7, 2.56e-7, 1.76e-7]

endf8_fispact = [2.75e-6, 8.7e-7, 6.79e-7, 6.12e-7, 4.63e-7]
eade = [2.95e-6, 8.92e-7, 7.18e-7, 6.18e-7, 4.68e-7]

endf8_icrp74 = np.array([2.19e-6, 8.20e-7, 5.80e-7, 5.21e-7, 3.73e-7])
experiment = np.array([2.46e-6, 6.99e-7, 4.95e-7, 4.16e-7, 3.16e-7])


fig, ax = plt.subplots()
ax.loglog(exp_data[:, 0], exp_data[:, 1], label='Experiment')
ax.loglog(times, endf8_icrp74, 'kx', label='OpenMC, ENDF/B-VIII.0, ICRP74')
#ax.loglog(times, eade, 'ro', label='Eade OpenMC (ICRP74)')
ax.set_xlabel('Time [d]')
ax.set_ylabel('Dose rate [Sv/h]')
ax.grid(which='both')
ax.legend()
plt.savefig('fng_dose_comparison.pdf', bbox_inches='tight')
#plt.show()


mcr2s_ce = [0.883, 1.115, 1.230, 1.233, 1.214]

fig, ax = plt.subplots()
plt.plot(range(5), endf8_icrp74/experiment, 'bx', label='OpenMC')
plt.plot(range(5), mcr2s_ce, 'rs', label='MCR2S')
xmin, xmax = -0.3, 4.3
plt.plot([xmin, xmax], [1., 1.], 'k--')
sigma = 0.1
verts = [(xmin, 1. - sigma), (xmin, 1. + sigma), (xmax, 1. + sigma), (xmax, 1. - sigma)]
poly = Polygon(verts, facecolor='gray', alpha=0.5)
ax.add_patch(poly)
ax.set_xlabel('Time', fontsize=14)
ax.set_ylabel('C/E', fontsize=14)
ax.set_ylim(0.7, 1.3)
ax.set_xticks(range(5), ['1 d', '7 d', '15 d', '30 d', '60 d'])
ax.set_xlim(xmin, xmax)
ax.legend()
ax.grid(which='both')
plt.savefig('fng_dose_ce.pdf', bbox_inches='tight')
#plt.show()
