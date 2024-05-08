import argparse
import json
from math import pi
from pathlib import Path
import warnings

import h5py
import openmc.deplete
import openmc
import dill
import numpy as np
from uncertainties import ufloat, UFloat


# Set cross sections
openmc.config['cross_sections'] = '/opt/data/hdf5/endfb-viii.0-hdf5/cross_sections.xml'

# Set chain for depletion / decay source generation
openmc.config['chain_file'] = 'chain_reduced.xml'

# Run arguments
nodes = 8
run_kwargs = {
    #'mpi_args': ['srun', '-N', f'{nodes}', '-n', f'{nodes*2}', '--cpu-bind=socket']
}


# TODO: Replace with activation_cells.json? Note that it currently includes
# cells 636, 637, 638, 648, 649, and 650 which are steel in the neutron model
# but turn into air in the photon model. Probably should be removed from
# activation_cells.json
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

# This doesn't include the source, but that is consistent with the cell-based
# workflow that we were doing before
ymin = 5.6

overall_bbox = openmc.BoundingBox(
    lower_left=(-49.5, 5.6, -49.2),
    upper_right=(49.5, 77.43, 49.2)
)
mesh_dimension = (20, 20, 20)


def fng_tabulated_source():
    # Load tabulated source data generated from compiled source
    with h5py.File('source_data.h5') as fh:
        energies = fh['energy_bins'][()]  # (3001,)
        thetas = fh['theta_bins'][()]     # (361,)
        yield_data = fh['yields'][()]     # (360, 3000)

    # Compute mu = cos(theta)
    cos_theta = np.cos(thetas)

    # Yield corresponding to each angle (summed over energies)
    intensity = yield_data.sum(axis=1)

    yields = []
    for i in range(intensity.size):
        yield_i = np.concatenate([yield_data[i], [0.0]])
        yields.append(openmc.stats.Tabular(energies, yield_i))

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

        # Create source for each angle
        source = openmc.IndependentSource(
            angle=angle_dist,
            energy=energy_dist,
            strength=intensity[i]
        )
        sources.append(source)

    return sources


def activation(model: openmc.Model, campaign: int, output_dir: Path):
    # Apply FNG source
    model.settings.source = fng_tabulated_source()

    if campaign == 1:
        hour = 3600.0
        day = 24*hour
        source_rates = [2.32e10, 0.0, 2.87e10, 0.0, 1.90e10, 0.0, 1.36e10]
        source_times = [19440., 61680., 32940., 54840., 15720., 6360., 8940.]
        cooling_times_cumulative = [
            1*hour, 6*hour, 12*hour, 16*hour, 20*hour, 1*day, 2*day, 3*day, 4*day,
            5*day, 7*day, 9*day, 12*day, 15*day, 18*day, 21*day, 30*day, 60*day
        ]
    elif campaign == 2:
        sources = [5.31e14, 3.35e14, 0., 9.50e14, 0., 1.29e14, 0., 4.0e12]
        source_times = [17480., 7820., 54140., 22140., 900., 3820., 420., 140.]
        source_rates = [s/t for s, t in zip(sources, source_times)]
        cooling_times_cumulative = [
            4380., 6180., 7488., 11580., 17280., 24480., 34080., 45780., 57240.,
            72550., 90720., 132000., 212400., 345600., 479300., 708500.,
            1050000., 1670000., 1710000.
        ]
    else:
        raise ValueError("Invalid FNG irradiation campaign")

    source_rates.extend([0.0]*len(cooling_times_cumulative))
    cooling_times = list(np.diff(cooling_times_cumulative, prepend=0.0))
    timesteps = source_times + cooling_times

    # Create mesh for activation / photon source
    mesh = openmc.RegularMesh()
    mesh.lower_left = overall_bbox.lower_left
    mesh.upper_right = overall_bbox.upper_right
    mesh.dimension = mesh_dimension

    # Get fluxes and micros based on cells
    fluxes, micros = openmc.deplete.get_microxs_and_flux(
        model, mesh, run_kwargs=run_kwargs)

    output_dir.mkdir(exist_ok=True)
    np.save(output_dir / 'fluxes.npy', fluxes)
    with open(output_dir / 'micros.pkl', 'wb') as fh:
        dill.dump(micros, fh)

    # Get homogenized materials
    activation_mats = mesh.get_homogenized_materials(model, include_void=False)

    # Save a copy of activation model
    model.export_to_model_xml(output_dir / 'activation_model.xml')

    # Create transport operator
    op = openmc.deplete.IndependentOperator(
        activation_mats, fluxes, micros, normalization_mode='source-rate')

    # Change output directory
    op.output_dir = output_dir

    # Run depletion
    predictor = openmc.deplete.PredictorIntegrator(op, timesteps, source_rates=source_rates)
    predictor.integrate(final_step=False)

    results = openmc.deplete.Results(op.output_dir / 'depletion_results.h5')

    # Get photon sources for each mesh element at each cooling time and serialize
    sources = {}
    for i_cool, _ in enumerate(cooling_times):
        sources[i_cool] = []
        for mat in activation_mats:
            mat = results[8+i_cool].get_material(str(mat.id))
            sources[i_cool].append(mat.decay_photon_energy)

    with open(op.output_dir / 'sources.pkl', 'wb') as fh:
        dill.dump(sources, fh)


def get_nickel_foil_rr(model: openmc.Model, output_dir: Path):
    # Add reaction rates in Ni58
    nickel_foil_tally = openmc.Tally()
    foil_cells = [602, 604, 608, 606, 605, 603]
    nickel_foil_tally.filters = [openmc.CellFilter(foil_cells)]
    nickel_foil_tally.nuclides = ['Ni58']
    nickel_foil_tally.scores = ['(n,2n)', '(n,p)']
    nickel_foil_tally.multiply_density = False

    model.tallies = [nickel_foil_tally]
    model.run(cwd=output_dir / 'nickel_foil_rr', **run_kwargs)


def photon_calculation(path_model: Path, campaign: int, dose_function: str, output_dir: Path):
    # Get Model object and add source and tallies
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", openmc.IDWarning)
        model = openmc.Model.from_model_xml(path_model)

    num_cooling_times = 18 if campaign == 1 else 19
    dose_values = []
    for cooling_time in range(num_cooling_times):
        generate_photon_sources(model, cooling_time, output_dir)
        generate_dose_tallies(model, dose_function)

        # Run OpenMC and compute dose
        intensity = sum(s.strength for s in model.settings.source)
        print(f'Photon transport calculation ({cooling_time}), source = {intensity:.3e} γ/s ...')
        sp_filename = model.run(output=False, cwd=output_dir / f'photon_{cooling_time}',
                                **run_kwargs)
        Sv_per_h = get_dose(sp_filename, campaign)
        print(f'Dose rate (flux) = {Sv_per_h} Sv/h')
        dose_values.append((Sv_per_h.nominal_value, Sv_per_h.std_dev))

    # Save dose values to numpy file
    np.save(output_dir / 'dose.npy', dose_values)


def get_dose(statepoint_file, campaign: int) -> UFloat:
    if campaign == 1:
        # Get volume of G-M dose meter (radius of 1.9 cm)
        r = 1.9
        cm3 = 4/3 * pi * r**3
    else:
        # Get volume of tissue-equivalent scintillator
        d = h = 4.6
        cm3 = pi * (d/2)**2 * h

    with openmc.StatePoint(statepoint_file) as sp:
        flux_tally = list(sp.tallies.values())[0]

    # Calculate dose using flux-to-dose conversion factor
    mean = flux_tally.mean.ravel()[0]
    stdev = flux_tally.std_dev.ravel()[0]
    pSv_cm3_per_sec = ufloat(mean, stdev)
    pSv_per_sec = pSv_cm3_per_sec / cm3
    Sv_per_h = pSv_per_sec * 1e-12 * 3600.
    return Sv_per_h


def generate_photon_sources(model, cooling_time: int, output_dir: Path):
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in dose_cell_ids]

    # Load photon energy distributions from activation calculation
    with open(output_dir / 'sources.pkl', 'rb') as fh:
        energy_dists = dill.load(fh)

    # Create source for each mesh element
    element_sources = []
    for energy in energy_dists[cooling_time]:
        source = openmc.IndependentSource(
            energy=energy,
            particle='photon',
            strength=energy.integral(),
        )
        element_sources.append(source)

    # Create mesh
    mesh = openmc.RegularMesh()
    mesh.lower_left = overall_bbox.lower_left
    mesh.upper_right = overall_bbox.upper_right
    mesh.dimension = mesh_dimension

    # Create mesh source
    model.settings.source = openmc.MeshSource(
        mesh, element_sources, constraints={'domains': dose_cells}
    )


def generate_dose_tallies(model: openmc.Model, dose_function: str):
    # Set dose function
    if dose_function == 'ans1977':
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
        # From Table A.21
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

    model.tallies = openmc.Tallies([flux_tally])


if __name__ == '__main__':
    # General configuration
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--campaign', type=int, choices=(1, 2), default=1)
    parser.add_argument('-n', '--model-neutron', type=Path, default='fng_neutron.xml')
    parser.add_argument('-p', '--model-photon', type=Path, default='fng_photon.xml')
    parser.add_argument('-f', '--dose-function', type=str, choices=('ans1977', 'icrp74', 'icrp116'), default='ans1977')
    parser.add_argument('-d', '--directory', type=Path, default=Path('results'))

    # Execution options
    parser.add_argument('--run-activation', action='store_true')
    parser.add_argument('--no-run-activation', action='store_false', dest='run_activation')
    parser.add_argument('--run-photon', action='store_true')
    parser.add_argument('--no-run-photon', action='store_false', dest='run_photon')
    parser.set_defaults(run_activation=True)
    parser.set_defaults(run_photon=True)
    args = parser.parse_args()

    # Step 1: Run neutron transport and activation
    if args.run_activation:
        neutron_model = openmc.Model.from_model_xml(args.model_neutron)
        activation(neutron_model, args.campaign, args.directory)

    # Step 2: Run photon transport for dose
    if args.run_photon:
        photon_calculation(args.model_photon, args.campaign, args.dose_function, args.directory)
