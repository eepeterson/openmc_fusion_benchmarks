import openmc
import numpy as np
import watts
import pypact as pp
import dill
import json

import fispact_utils


with open('activation_cells.json', 'r') as fh:
    dose_cell_ids = json.load(fh)
with open('bounding_boxes.json', 'r') as fh:
    bounding_boxes = json.load(fh)

openmc.config['cross_sections'] = '/opt/data/hdf5/endfb80_hdf5/cross_sections.xml'


def calculate_volumes(cell_ids):
    model = openmc.Model.from_xml()
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in dose_cell_ids]

    # Run a volume calculation
    vol_calcs = []
    for cell in dose_cells:
        lower_left, upper_right = bounding_boxes[str(cell.id)]
        if np.isinf(lower_left).any() or np.isinf(upper_right).any():
            raise ValueError(f'Cell {cell.id} has an infinite bounding box')
        vol_calcs.append(openmc.VolumeCalculation([cell], 1_000_000, lower_left, upper_right))

    model.settings.volume_calculations = vol_calcs
    model.calculate_volumes(threads=16)

    volumes = {cell.id: cell.volume for cell in dose_cells}
    with open('cell_volumes.json', 'w') as fh:
        json.dump(volumes, fh)


def apply_volumes(filename):
    # Read volumes from JSON
    with open(filename, 'r') as fh:
        volumes = json.load(fh)

    # Add volumes to cells
    model = openmc.Model.from_xml()
    cells = model.geometry.get_all_cells()
    for uid, volume in volumes.items():
        cells[int(uid)].volume = volume

    # Re-export model
    model.export_to_xml()


def neutron_flux(dose_cell_ids):
    model = openmc.Model.from_xml()

    tally = openmc.Tally()
    tally.filters = [
        openmc.CellFilter(dose_cell_ids),
        openmc.EnergyFilter.from_group_structure('CCFE-709'),
    ]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    # Run OpenMC
    sp_path = model.run(threads=12)

    # Get flux for one specific material
    with openmc.StatePoint(sp_path) as sp:
        return sp.tallies[tally.id].get_reshaped_data().squeeze()


def fispact(dose_cell_ids, flux):
    day = 24*3600
    source_rate_irradiate = [2.32e10, 0.0, 2.87e10, 0.0, 1.90e10, 0.0, 1.36e10]
    timesteps_irradiate = [19440., 61680., 32940., 54840., 15720., 6360., 8940.]
    cooling_times = [1*day, 6*day, 8*day, 15*day, 30*day]

    model = openmc.Model.from_xml()
    cells = model.geometry.get_all_cells()

    sources = {cooling_time: {} for cooling_time in [1, 7, 15, 30, 60]}

    for i, target_cell_id in enumerate(dose_cell_ids):
        neutron_cm_per_src = flux[i]
        src_per_sec = source_rate_irradiate[0]

        material = cells[target_cell_id].fill
        cm3 = material.volume
        neutron_per_cm2_sec = neutron_cm_per_src * src_per_sec / cm3

        # FISPACT-II requires that fluxes be in descending energy
        # See https://fispact.ukaea.uk/forum/viewtopic.php?t=30
        fluxes = np.concatenate((neutron_per_cm2_sec[::-1].copy(), [1.0]))
        np.savetxt('fluxes', fluxes, footer=f'OpenMC flux in cell {target_cell_id}', comments='')

        # Create FISPACT watts plugin
        fispact = watts.PluginGeneric(
            executable='fispact',
            execute_command='{self.executable} {self.base_name}',
            template_file='inventory.tmpl',
            extra_inputs=['files', 'fluxes', 'COLLAPX', 'ARRAYX'],
            show_stdout=True, show_stderr=True
        )
        fispact.input_name = 'openmc_compare.i'
        fispact.base_name = fispact.input_name[:-2]

        atoms = material.get_nuclide_atoms()
        fluxes = neutron_cm_per_src.sum() / cm3 * np.asarray(source_rate_irradiate)

        params = watts.Parameters()
        params['density'] = material.get_mass_density()
        params['num_nuclides'] = len(atoms)
        params['atoms'] = '\n'.join(f'{nuc} {atom:.6e}' for nuc, atom in atoms.items())
        params['irradiation_times'] = watts.Quantity(timesteps_irradiate, 's')
        params['fluxes'] = fluxes
        params['cooling_times'] = watts.Quantity(cooling_times, 's')
        params['cell_id'] = target_cell_id

        result = fispact(params, name=f'FISPACT, cell {target_cell_id}')

        # Get gamma spectrum and write to source
        with pp.Reader(result.base_path / f'openmc_compare.out') as output:
            ...

        for i_cool, cooling_time in enumerate([1, 7, 15, 30, 60]):
            # Generate distribution
            timestep = output[-5 + i_cool]
            sources[cooling_time][target_cell_id] = fispact_utils.fispact_photon_energy(
                timestep, material.get_mass_density()
            )

    with open('sources_fispact.pkl', 'wb') as fh:
        dill.dump(sources, fh)


#calculate_volumes(dose_cell_ids)
#apply_volumes('cell_volumes.json')
flux = neutron_flux(dose_cell_ids)
