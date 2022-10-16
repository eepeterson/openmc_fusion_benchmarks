import argparse
from math import pi

import openmc
import openmc.data
import numpy as np
import dill
from uncertainties import ufloat

from dose_cells import dose_cell_ids, inner_cell_ids, next_cell_ids, front_cell_ids


def main(cooling_time):
    generate_sources(cooling_time)
    generate_tallies()
    print('Running OpenMC...')
    openmc.run(output=False)
    Sv_per_h = get_dose('statepoint.40.h5')
    print(f'Dose rate (flux) = {Sv_per_h} Sv/h')


def get_dose(statepoint_file):
    # Get volume of cell 651 (radius of 1.9 cm)
    r = 1.9
    cm3 = 4/3 * pi * r**3
    g_per_cm3 = 0.0011050927231898923
    kg = g_per_cm3*cm3 * 1e-3

    with openmc.StatePoint(statepoint_file) as sp:
        tallies = list(sp.tallies.values())
        flux_tally, heating_tally = tallies

    # Calculate dose using flux-to-dose conversion factor
    mean = flux_tally.mean.ravel()[0]
    stdev = flux_tally.std_dev.ravel()[0]
    mrem_cm3_per_sec = ufloat(mean, stdev)
    mrem_per_sec = mrem_cm3_per_sec / cm3
    Sv_per_mrem = 1e-5
    Sv_per_h = mrem_per_sec * 3600. * Sv_per_mrem
    return Sv_per_h

    """
    # Calculate dose using energy deposition directly
    mean = heating_tally.mean.ravel()[0]
    stdev = heating_tally.std_dev.ravel()[0]
    eV_per_sec = ufloat(mean, stdev)
    J_per_eV = 1.602176634e-19
    J_per_kg_per_sec = eV_per_sec / kg * J_per_eV
    Sv_per_h = J_per_kg_per_sec * 3600.
    print(f'Dose rate (Edep) = {Sv_per_h} Sv/h')
    """


def generate_sources(cooling_time: int, cutoff=1100.0):
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
        energy = energy_dists[cooling_time][cell.id]

        # Restrict energies to > 1.1 keV
        mask = energy.x > cutoff
        energy.x = energy.x[mask]
        energy.p = energy.p[mask]

        source = openmc.Source(
            space=space,
            energy=energy,
            particle='photon',
            strength=energy.integral(),
            domains=[cell]
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


def generate_tallies():
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

    tallies = openmc.Tallies([flux_tally, heating_tally])
    tallies.export_to_xml()


if __name__ == '__main__':
    cooling_times = [1, 7, 15, 30, 60]
    parser = argparse.ArgumentParser()
    parser.add_argument('cooling_time', type=int, choices=cooling_times)
    args = parser.parse_args()
    main(args.cooling_time)
