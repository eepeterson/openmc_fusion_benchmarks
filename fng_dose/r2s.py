import argparse
import json
from math import pi
from pathlib import Path
import warnings

import openmc.deplete
import openmc
import dill
import numpy as np
from uncertainties import ufloat


# Set cross sections
openmc.config['cross_sections'] = '/opt/data/hdf5/endfb-viii.0-hdf5/cross_sections.xml'

# Set chain for depletion / decay source generation
openmc.config['chain_file'] = 'chain_reduced.xml'

# TODO: Replace with activation_cells.json?
inner_cell_ids = [
    160, 161, 173, 174, 175, 176, 177, 178, 179, 226, 228, 230, 239, 240, 242,
    251, 252, 253, 601, 620, 621, 622, 624, 625, 627, 628
]
next_cell_ids = [
    137, 138, 139, 140, 162, 163, 180, 181, 231, 232, 243, 244, 254, 255, 256,
    274, 275, 276, 277, 278
]
front_cell_ids = [
    102, 103, 104, 105, 106, 107, 119, 142, 165, 114, 115, 116, 117, 118, 141,
    164, 182, 233, 245
]
dose_cell_ids = inner_cell_ids + next_cell_ids + front_cell_ids


def calculate_volumes(cell_ids):
    model = openmc.Model.from_xml()
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in cell_ids]

    # Get bounding boxes
    with open('bounding_boxes.json', 'r') as fh:
        bounding_boxes = json.load(fh)

    # Run a volume calculation
    vol_calcs = []
    for cell in dose_cells:
        lower_left, upper_right = bounding_boxes[str(cell.id)]
        if np.isinf(lower_left).any() or np.isinf(upper_right).any():
            raise ValueError(f'Cell {cell.id} has an infinite bounding box')
        vol_calcs.append(openmc.VolumeCalculation([cell], 1_000_000, lower_left, upper_right))

    model.settings.volume_calculations = vol_calcs
    model.calculate_volumes()

    volumes = {cell.id: cell.volume for cell in dose_cells}
    with open('cell_volumes.json', 'w') as fh:
        json.dump(volumes, fh)


def apply_volumes(model, filename='cell_volumes.json', material=False):
    # Read volumes from JSON
    with open(filename, 'r') as fh:
        volumes = json.load(fh)

    # Add volumes to cells/materials
    cells = model.geometry.get_all_cells()
    for uid, volume in volumes.items():
        cell = cells[int(uid)]
        if material:
            cell.fill = cell.fill.clone()
            cell.fill.depletable = True
            cell.fill.volume = volume
        else:
            cell.volume = volume

    # Make sure new materials in geometry get exported
    if material:
        model.materials.clear()


def activation(model: Path, operator_type: str, output_dir: Path):
    model = openmc.Model.from_model_xml(model)

    schedule = 'campaign1'

    if schedule == 'campaign1':
        day = 24*3600
        source_rates = [
            2.32e10, 0.0, 2.87e10, 0.0, 1.90e10, 0.0, 1.36e10,
            0.0, 0.0, 0.0, 0.0, 0.0
        ]
        timesteps = [
            19440., 61680., 32940., 54840., 15720., 6360., 8940.,
            1*day, 6*day, 8*day, 15*day, 30*day
        ]
    elif schedule == 'campaign2':
        sources = [5.31e14, 3.35e14, 0., 9.50e14, 0., 1.29e14, 0., 4.0e12]
        timesteps = [17480., 7820., 54140., 22140., 900., 3820., 420., 140.]
        source_rates = np.array(sources) / np.array(timesteps)
    elif schedule == 'eff726':
        # Cooling time of 1 min
        source_rates = [1.0e11, 0.0]
        timesteps = [1e4, 60.0]
    else:
        raise ValueError("Invalid schedule")

    # Get list of cells for activation
    all_cells = model.geometry.get_all_cells()
    cells = [all_cells[uid] for uid in dose_cell_ids]

    if operator_type == 'independent':
        # Apply volumes to cells
        apply_volumes(model, material=False)

        # Get fluxes and micros based on cells
        fluxes, micros = openmc.deplete.get_microxs_and_flux(model, cells)

        output_dir.mkdir(exist_ok=True)
        np.save(output_dir / 'fluxes.npy', fluxes)
        with open(output_dir / 'micros.pkl', 'wb') as fh:
            dill.dump(micros, fh)

        # Create material copies for activation
        activation_mats = openmc.Materials()
        activation_mat_by_id = {}
        for cell in cells:
            mat = cell.fill.clone()
            mat.name = f'Cell {cell.id}'
            mat.depletable = True
            mat.volume = cell.volume
            activation_mats.append(mat)
            activation_mat_by_id[cell.id] = mat

        # Create transport operator
        op = openmc.deplete.IndependentOperator(
            activation_mats, fluxes, micros, normalization_mode='source-rate')

    else:
        # Apply volumes to materials
        apply_volumes(model, material=True)

        # Get dictionary of materials
        activation_mat_by_id = {c.id: c.fill for c in cells}

        # Create transport operator
        op = openmc.deplete.CoupledOperator(model, normalization_mode='source-rate')

    # Change output directory
    op.output_dir = output_dir

    # Run depletion
    predictor = openmc.deplete.PredictorIntegrator(op, timesteps, source_rates=source_rates)
    predictor.integrate(final_step=False)

    sources = {}
    results = openmc.deplete.Results(op.output_dir / 'depletion_results.h5')
    for i_cool, cooling_time in enumerate([1, 7, 15, 30, 60]):
        sources[cooling_time] = {}
        for uid in dose_cell_ids:
            mat_id = activation_mat_by_id[all_cells[uid].id].id
            mat = results[8+i_cool].get_material(str(mat_id))
            sources[cooling_time][uid] = mat.decay_photon_energy

    with open(op.output_dir / 'sources.pkl', 'wb') as fh:
        dill.dump(sources, fh)


def photon_calculation(path_model: Path, cooling_time, dose_function: str, output_dir: Path):
    # Get Model object and add source and tallies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", openmc.IDWarning)
        model = openmc.Model.from_model_xml(path_model)
    generate_photon_sources(model, cooling_time, output_dir)
    generate_dose_tallies(model, dose_function)

    # Run OpenMC and compute dose
    intensity = sum(s.strength for s in model.settings.source)
    print(f'Photon transport calculation at {cooling_time} d, source = {intensity:.3e} γ/s ...')
    sp_filename = model.run(output=False, cwd=output_dir / f'photon_{cooling_time}d')
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


def generate_photon_sources(model, cooling_time: int, output_dir: Path):
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in dose_cell_ids]

    # Load photon energy distributions from activation calculation
    with open(output_dir / 'sources.pkl', 'rb') as fh:
        energy_dists = dill.load(fh)

    # Load bounding box information
    with open('bounding_boxes.json', 'r') as fh:
        bounding_boxes = json.load(fh)

    # Create sources for each depletable region
    sources = []
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

    model.settings.source = sources


def generate_dose_tallies(model: openmc.Model, dose_function: str):
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

    # General configuration
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--model-neutron', type=Path, default='fng_neutron.xml')
    parser.add_argument('-p', '--model-photon', type=Path, default='fng_photon.xml')
    parser.add_argument('-o', '--operator', choices=('coupled', 'independent'), default='independent')
    parser.add_argument('-f', '--dose-function', type=str, choices=('ans1977', 'icrp74', 'icrp116', 'ethan'), default='ans1977')
    parser.add_argument('-d', '--directory', type=Path, default=Path('results'))
    parser.add_argument('cooling_times', type=int, nargs='*', default=cooling_times)

    # Execution options
    parser.add_argument('--run-volume-calc', action='store_true')
    parser.add_argument('--no-run-volume-calc', action='store_false', dest='run_volume_calc')
    parser.add_argument('--run-activation', action='store_true')
    parser.add_argument('--no-run-activation', action='store_false', dest='run_activation')
    parser.add_argument('--run-photon', action='store_true')
    parser.add_argument('--no-run-photon', action='store_false', dest='run_photon')
    parser.set_defaults(run_volume_calc=False)
    parser.set_defaults(run_activation=True)
    parser.set_defaults(run_photon=True)
    args = parser.parse_args()

    # Preprocessing step: calculate volumes
    if args.run_volume_calc:
        with open('activation_cells.json', 'r') as fh:
            dose_cell_ids = json.load(fh)
        calculate_volumes(dose_cell_ids)

    # Step 1: Run neutron transport and activation
    if args.run_activation:
        activation(args.model_neutron, args.operator, args.directory)

    # Step 2: Run photon transport for dose
    if args.run_photon:
        for time in args.cooling_times:
            if time not in cooling_times:
                raise ValueError(f"Invalid cooling time: {time}")
            photon_calculation(args.model_photon, time, args.dose_function, args.directory)
