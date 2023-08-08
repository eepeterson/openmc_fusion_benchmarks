import numpy as np
import json
import h5py

import openmc

sp = openmc.StatePoint('statepoint.50.h5')
flux_tally = sp.get_tally(name='flux_tally')
cell_ids = flux_tally.find_filter(openmc.CellFilter).bins
ebins = flux_tally.find_filter(openmc.EnergyFilter).values
flux = flux_tally.get_reshaped_data().squeeze()
print(flux.shape)
rel_err = flux_tally.get_reshaped_data(value='rel_err').squeeze()
cells = sp.summary.geometry.get_all_cells()

with open('../cell_volumes.json', 'r') as fh:
    volume_dict = json.load(fh)
volumes = np.array([volume_dict[str(uid)] for uid in cell_ids])
print(volumes)
flux /= volumes[:, np.newaxis]

flux_map = {}
flux_error_map = {}
for i, uid in enumerate(cell_ids):
    flux_map[uid] = flux[i, :]
    flux_error_map[uid] = rel_err[i, :]

with h5py.File('flux_spectrum_results.h5', 'a') as fh:
    openmc_grp = fh.create_group('openmc')
    openmc_grp.create_dataset('ebin_edges', data=ebins)
    for cell_id in flux_map.keys():
        tmp_grp = openmc_grp.create_group(f"cell_{cell_id}")
        tmp_grp.create_dataset('flux', data=flux_map[cell_id])
        tmp_grp.create_dataset('rel_err', data=flux_error_map[cell_id])
