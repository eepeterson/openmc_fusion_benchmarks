import openmc
import openmc.data
import numpy as np
import dill

from dose_cells import dose_cell_ids, inner_cell_ids, next_cell_ids, front_cell_ids

model = openmc.Model.from_xml()
cells = model.geometry.get_all_cells()
dose_cells = [cells[uid] for uid in dose_cell_ids]

with open('../neutron/sources.pkl', 'rb') as fh:
    energy_dists = dill.load(fh)

# Create sources for each depletable region
sources = []
intensity_inner = 0.0
intensity_next = 0.0
intensity_front = 0.0
for cell in dose_cells:
    space = openmc.stats.Box(*cell.bounding_box)
    energy = energy_dists[cell.id]

    source = openmc.Source(
        space=space,
        energy=energy,
        particle='photon',
        strength=energy.integral(),
        cells=[cell]
    )
    sources.append(source)

    if cell.id in inner_cell_ids:
        intensity_inner += source.strength
    elif cell.id in next_cell_ids:
        intensity_next += source.strength
    elif cell.id in front_cell_ids:
        intensity_front += source.strength
model.settings.source = sources
model.settings.export_to_xml()

print(f'Source (inner) = {intensity_inner} γ/s')
print(f'Source (next)  = {intensity_next} γ/s')
print(f'Source (front) = {intensity_front} γ/s')

# Dose based on ICRU flux-to-dose
energies = 1e6*np.array(
    [0.001, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3,
     0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
mrem_per_hour = [7.43e-07, 3.12e-07,1.68e-07,7.21e-08,4.29e-08,3.23e-08,2.89e-08,3.07e-08,
                 3.71e-08,5.99e-08,8.56e-08, 1.38e-07,1.893e-07,2.38e-07,2.84e-07,3.69e-07,4.47e-07,
                 6.14e-07,7.55e-07,9.96e-07,1.21e-06,1.41e-06,1.61e-06,2.01e-06,2.4e-06,2.4e-06]
dose_func = energies, mrem_per_hour
dose_filter = openmc.EnergyFunctionFilter(*dose_func)
particle_filter = openmc.ParticleFilter(['photon'])
cell_filter = openmc.CellFilter([651])
flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter, dose_filter, particle_filter]
flux_tally.scores = ['flux']

# Dose based on energy deposition
heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter]
heating_tally.scores = ['heating']

model.tallies = openmc.Tallies([flux_tally, heating_tally])
model.tallies.export_to_xml()
