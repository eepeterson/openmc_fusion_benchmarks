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
print(elem_width)

plot = openmc.lib.plot._PlotBase()
plot.width = width[0]
plot.height = width[1]
plot.h_res = res
plot.v_res = res
plot.basis = 'xy'

#with open('activation_cells.json', 'r') as fh:
with open('valid_cells.json', 'r') as fh:
    dose_cell_ids = json.load(fh)


def get_clean_ll_ur(lower_left, upper_right, cell):
    ll_ref, ur_ref = cell.bounding_box
    print(f"cell {cell.id} auto bb: {ll_ref, ur_ref}")
    if np.isinf(ll_ref).any() or np.isinf(ur_ref).any():
        return lower_left, upper_right
    if (np.any(np.asarray(lower_left) > ll_ref) or
        np.any(np.asarray(upper_right) < ur_ref)):
        return ll_ref, ur_ref
    else:
        return lower_left, upper_right


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

    # Expand bounding box and make sure no width is 0
    new_width = new_ur - new_ll + 2*elem_width
    new_ll -= expand*new_width
    new_ur += expand*new_width

    return new_ll, new_ur

# Get 3D array of cell IDs in (x, y, z) order 
z_array = (np.linspace(lower_left[2], upper_right[2], res + 1) + 0.5*elem_width[2])[:-1]
cell_ids = []
for z in tqdm(z_array):
    plot.origin = (center[0], center[1], z)
    img_data = openmc.lib.id_map(plot)
    cell_ids_z = np.swapaxes(img_data[..., 0], 0, 1)
    cell_ids.append(cell_ids_z[:, ::-1])
cell_ids = np.dstack(cell_ids)

bounding_boxes = {}
for uid in dose_cell_ids:
    # Need to check a few cases because of bug and small volumes
    if uid == 4:
        ll, ur = [-3.0, 0.0, -3.0], [3.0, 0.3, 3.0]
    elif uid == 46:
        ll, ur = [-17.7, -1.3, -17.7], [17.7, 5.6, 17.7]
    elif uid == 47:
        ll, ur = [-17.7, -1.3, -17.7], [17.7, 5.6, 17.7]
    elif uid == 48:
        ll, ur = [-17.7, -13.5, -17.7], [17.7, -1.3, 17.7]
    elif uid == 49:
        ll, ur = [-17.7, -13.5, -17.7], [17.7, -1.3, 17.7]
    elif uid == 608:
        ll, ur = [7.1, 33.0, -1.0], [7.4, 34.0, 1.0]
    else:
        ll, ur = get_new_bb2(cell_ids, lower_left, upper_right, res, uid,
                             expand=1.)

    ll, ur = get_clean_ll_ur(ll, ur, openmc.lib.cells[uid])
    print(f"cell {uid}:", ll, ur)
    bounding_boxes[uid] = (tuple(ll), tuple(ur))

with open('bounding_boxes.json', 'w') as fh:
    json.dump(bounding_boxes, fh)

openmc.lib.finalize()
