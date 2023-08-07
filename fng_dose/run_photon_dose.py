import argparse
import json
from math import pi
from pathlib import Path

import openmc
import openmc.data
import numpy as np
import dill
from uncertainties import ufloat

from dose_cells import dose_cell_ids, inner_cell_ids, next_cell_ids, front_cell_ids
import data_config


def main(path_model, cooling_time, dose_function: str):
    model = openmc.Model.from_model_xml(path_model)
    generate_sources(model, cooling_time)
    generate_tallies(model, dose_function)
    print('Running OpenMC...')
    sp_filename = model.run(output=False, cwd='photon_dose')
    Sv_per_h = get_dose(sp_filename)
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
    pSv_cm3_per_sec = ufloat(mean, stdev)
    pSv_per_sec = pSv_cm3_per_sec / cm3
    Sv_per_h = pSv_per_sec * 1e-12 * 3600.
    return Sv_per_h

    """
    # Calculate dose using energy deposition directly
    mean = heating_tally.mean.ravel()[0]
    stdev = heating_tally.std_dev.ravel()[0]
    eV_per_sec = ufloat(mean, stdev)
    J_per_eV = 1.602176634e-19
    J_per_kg_per_sec = eV_per_sec / kg * J_per_eV
    Sv_per_h_energy = J_per_kg_per_sec * 3600.
    print(f'Dose rate (Edep) = {Sv_per_h_energy} Sv/h')
    """


def generate_sources(model, cooling_time: int):
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in dose_cell_ids]

    # Load photon energy distributions from activation calculation
    with open('activation/sources.pkl', 'rb') as fh:
        energy_dists = dill.load(fh)

    # Load bounding box information
    with open('bounding_boxes.json', 'r') as fh:
        bounding_boxes = json.load(fh)


    # Create sources for each depletable region
    sources = []
    intensity_inner = 0.0
    intensity_next = 0.0
    intensity_front = 0.0
    for cell in dose_cells:
        space = openmc.stats.Box(*bounding_boxes[str(cell.id)])
        energy = energy_dists[cooling_time][cell.id]

        source = openmc.IndependentSource(
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

    print(f'Source (inner) = {intensity_inner} γ/s')
    print(f'Source (next)  = {intensity_next} γ/s')
    print(f'Source (front) = {intensity_front} γ/s')


def generate_tallies(model, dose_function: str):
    # Set dose function
    if dose_function == 'ethan':
        energies = 1e6*np.array([
            0.001, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3,
            0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0
        ])
        # Values are mrem-cm²
        pSv_cm2 = 1e12 * 1e-5 * np.array([
            7.43e-07, 3.12e-07,1.68e-07,7.21e-08,4.29e-08,3.23e-08,2.89e-08,3.07e-08,
            3.71e-08,5.99e-08,8.56e-08, 1.38e-07,1.893e-07,2.38e-07,2.84e-07,3.69e-07,4.47e-07,
            6.14e-07,7.55e-07,9.96e-07,1.21e-06,1.41e-06,1.61e-06,2.01e-06,2.4e-06,2.4e-06
        ])
    elif dose_function == 'ans1977':
        energies = 1e6*np.array([
            0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
            0.55, 0.6, 0.65, 0.7, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 2.8, 3.25, 3.75, 4.25,
            4.75, 5.0, 5.25, 5.75, 6.25, 6.75, 7.5, 9.0, 11.0, 13.0, 15.0
        ])
        # Values are rem-cm²-sec/hr
        pSv_cm2 = 1e-2 * 1e12 / 3600 * np.array([
            3.96e-6, 5.82e-7, 2.9e-7, 2.58e-7, 2.83e-7, 3.79e-7, 5.01e-7, 6.31e-7,
            7.59e-7, 8.78e-7, 9.85e-7, 1.08e-6, 1.17e-6, 1.27e-6, 1.36e-6, 1.44e-6,
            1.52e-6, 1.68e-6, 1.98e-6, 2.51e-6, 2.99e-6, 3.42e-6, 3.82e-6, 4.01e-6,
            4.41e-6, 4.83e-6, 5.23e-6, 5.6e-6, 5.8e-6, 6.01e-6, 6.37e-6, 6.74e-6,
            7.11e-6, 7.66e-6, 8.77e-6, 1.03e-5, 1.18e-5, 1.33e-5
        ])
    elif dose_function == 'icrp74':
        energies = 1e6*np.array([
            0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2,
            0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0
        ])
        pSv_cm2 = np.array([
            0.061, 0.83, 1.05, 0.81, 0.64, 0.55, 0.51, 0.53, 0.61, 0.89, 1.2, 1.8, 2.38,
            2.93, 3.44, 4.38, 5.20, 6.9, 8.6, 11.1, 13.4, 15.5, 17.6, 21.6, 25.6
        ])
    elif dose_function == 'icrp116':
        energies, pSv_cm2 = openmc.data.dose_coefficients('photon')

    dose_filter = openmc.EnergyFunctionFilter(energies, pSv_cm2)
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


if __name__ == '__main__':
    cooling_times = [1, 7, 15, 30, 60]
    parser = argparse.ArgumentParser()
    parser.add_argument('cooling_time', type=int, choices=cooling_times)
    parser.add_argument('-m', '--model', type=Path, default='fng_photon.xml')
    parser.add_argument('-f', '--function', type=str, choices=('ans1977', 'icrp74', 'icrp116', 'ethan'), default='ans1977')
    args = parser.parse_args()
    main(args.model, args.cooling_time, args.function)
