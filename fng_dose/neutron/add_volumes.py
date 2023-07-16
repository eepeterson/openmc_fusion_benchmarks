import json

import openmc
import numpy as np

import data_config


with open('activation_cells.json', 'r') as fh:
    dose_cell_ids = json.load(fh)
with open('bounding_boxes.json', 'r') as fh:
    bounding_boxes = json.load(fh)


def calculate_volumes(cell_ids):
    model = openmc.Model.from_xml()
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in cell_ids]

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

    # Re-export model
    model.export_to_xml()


if __name__ == '__main__':
    calculate_volumes(dose_cell_ids)
