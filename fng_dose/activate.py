import argparse
from pathlib import Path

import openmc.deplete
import openmc
import dill
import numpy as np

from dose_cells import dose_cell_ids
from volumes import apply_volumes
import data_config


def activation(model: Path, operator_type: str):
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
    op.output_dir = 'activation'

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

    with open('activation/sources.pkl', 'wb') as fh:
        dill.dump(sources, fh)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--model-neutron', type=Path, default='fng_neutron.xml')
    parser.add_argument('-o', '--operator', choices=('coupled', 'independent'), default='independent')
    parser.add_argument('--activation', action='store_true')
    parser.add_argument('--no-activation', dest='activation', action='store_false')
    parser.set_defaults(activation=True)
    args = parser.parse_args()

    # Run activation
    if args.activation:
        activation(args.model_neutron, args.operator)
