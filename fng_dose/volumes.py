import json

import openmc
import numpy as np

#import data_config


with open('valid_cells.json', 'r') as fh:
    dose_cell_ids = json.load(fh)

with open('bounding_boxes.json', 'r') as fh:
    bounding_boxes = json.load(fh)

print(len(bounding_boxes))


def calculate_volumes(cell_ids):
    model = openmc.Model.from_model_xml()
    cells = model.geometry.get_all_cells()
    dose_cells = [cells[uid] for uid in cell_ids]

    # Run a volume calculation
    vol_calcs = []
    for cell in dose_cells:
        lower_left, upper_right = bounding_boxes[str(cell.id)]
        vcalc = openmc.VolumeCalculation([cell], 1_000_000, lower_left, upper_right)
        vcalc.set_trigger(1e-3, 'rel_err')
        vol_calcs.append(vcalc)

    model.settings.volume_calculations = vol_calcs
    model.export_to_model_xml()
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


if __name__ == '__main__':
    calculate_volumes(dose_cell_ids)
