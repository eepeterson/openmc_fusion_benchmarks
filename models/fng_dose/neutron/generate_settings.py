from math import pi

import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.stats

# Take discrete values and create intervals that cover [-1, 1]
cos_theta = np.array([
    1.00000, 0.98481, 0.93969, 0.86603, 0.76604, 0.64279, 0.50000, 0.34202,
    0.17365, 0.00000, -0.17365, -0.34202, -0.50000, -0.64279, -0.76604,
    -0.86603, -0.93969, -0.98481, -1.00000
])
cos_theta = np.concatenate([[1.], cos_theta[:-1] + 0.5*np.diff(cos_theta), [-1.]])

# Weight relative intensity by width of corresponding cos_theta intervals
relative_intensity = np.array([
  1.05459, 1.05374, 1.05124, 1.04715, 1.04162, 1.03484, 1.02701, 1.01842,
  1.00932, 1.00000, 0.99075, 0.98184, 0.97354, 0.96609, 0.95969, 0.95452,
  0.95073, 0.94842, 0.94764
])
relative_intensity *= np.diff(-cos_theta)

yield_data = np.loadtxt('yields.txt')
yields = []
for i in range(relative_intensity.size):
    yields.append(
        openmc.stats.Tabular(1e6*yield_data[:, 2*i], yield_data[:, 2*i+1])
    )

sources = []
for i, (mu_high, mu_low) in enumerate(zip(cos_theta[:-1], cos_theta[1:])):
    mu_dist = openmc.stats.Uniform(mu_low, mu_high)
    phi_dist = openmc.stats.Uniform(0., 2*pi)
    angle_dist = openmc.stats.PolarAzimuthal(
        mu_dist,
        phi_dist,
        reference_uvw=(0.0, 1.0, 0.0),
    )
    energy_dist = yields[i]

    # Create sources so that strengths add up to 1
    source = openmc.Source(
        angle=angle_dist,
        energy=energy_dist,
        strength=relative_intensity[i] / relative_intensity.sum()
    )
    sources.append(source)

settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.particles = 10_000
settings.batches = 100
settings.source = sources
settings.export_to_xml()

"""
fig, ax = plt.subplots()
for idx, label in zip([0, 6, 9, 12, 18], [0, 60, 90, 120, 180]):
    E = yields[idx].x
    y = yields[idx].p
    ax.plot(E[:-1], y[:-1] / np.abs(np.diff(E)), label=f'{label} degree')
ax.set_xlabel('Energy [MeV]')
ax.set_ylabel('Neutron spectrum [n/MeV]')
ax.legend()
plt.show()
"""
