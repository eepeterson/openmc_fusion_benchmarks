import json

import openmc.lib
import numpy as np
from tqdm import tqdm

openmc.lib.init()

lower_left = np.array([-60.0, -15., -55.])
upper_right = np.array([60.0, 80., 55.])
center = (lower_left + upper_right) / 2
width = upper_right - lower_left
res = 500
elem_width = (upper_right - lower_left) / res

plot = openmc.lib.plot._PlotBase()
plot.width = width[0]
plot.height = width[1]
plot.h_res = res
plot.v_res = res
plot.basis = 'xy'

#with open('activation_cells.json', 'r') as fh:
with open('valid_cells.json', 'r') as fh:
    dose_cell_ids = json.load(fh)


def get_new_bb2(cell_ids, lower_left, upper_right, res, target_cell_id, expand=0.2):
    # Construct x, y, and z grids for this view
    elem_width = (upper_right - lower_left) / res
    x_grid = (np.linspace(lower_left[0], upper_right[0], res + 1) + 0.5*elem_width[0])[:-1]
    y_grid = (np.linspace(lower_left[1], upper_right[1], res + 1) + 0.5*elem_width[1])[:-1]
    z_grid = (np.linspace(lower_left[2], upper_right[2], res + 1) + 0.5*elem_width[2])[:-1]

    # Get set of points (N, 3) where this cell was found
    ix, iy, iz = np.where(cell_ids == target_cell_id)
    xs = x_grid[ix]
    ys = y_grid[iy]
    zs = z_grid[iz]
    points = np.column_stack((xs, ys, zs))

    new_ll = np.min(points, axis=0)
    new_ur = np.max(points, axis=0)

    # Expand bounding box
    new_width = new_ur - new_ll
    new_ll -= expand*new_width
    new_ur += expand*new_width

    return new_ll, new_ur

# Get 3D array of cell IDs in (x, y, z) order
z_array = (np.linspace(lower_left[2], upper_right[2], res + 1) + elem_width[2])[:-1]
cell_ids = []
for z in tqdm(z_array):
    plot.origin = (center[0], center[1], z)
    img_data = openmc.lib.id_map(plot)
    cell_ids_z = np.swapaxes(img_data[..., 0], 0, 1)
    cell_ids.append(cell_ids_z[:, ::-1])
cell_ids = np.dstack(cell_ids)

bounding_boxes = {}
for uid in dose_cell_ids:
    if not (cell_ids == uid).any():
        print(f'No match for cell {uid}')
        continue
    ll, ur = get_new_bb2(cell_ids, lower_left, upper_right, res, uid)
    print(uid, ll, ur)
    bounding_boxes[uid] = (tuple(ll), tuple(ur))

with open('bounding_boxes2.json', 'w') as fh:
    json.dump(bounding_boxes, fh)

openmc.lib.finalize()
