import openmc.deplete
import openmc
import dill
import numpy as np

from dose_cells import dose_cell_ids

model = openmc.Model.from_xml()

schedule = 'campaign1'

if schedule == 'campaign1':
    source_rates = [2.32e10, 0.0, 2.87e10, 0.0, 1.90e10, 0.0, 1.36e10, 0.0]
    timesteps = [19440., 61680., 32940., 54840., 15720., 6360., 8940., 7*24*3600.]
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

model.deplete(
    timesteps,
    source_rates=source_rates,
    method='predictor',
    final_step=False,
    operator_kwargs={
        'chain_file': 'chain_reduced.xml',
        'normalization_mode': 'source-rate',
    }
)

sources = {}
cells = model.geometry.get_all_cells()
results = openmc.deplete.Results('depletion_results.h5')
for uid in dose_cell_ids:
    mat_id = cells[uid].fill.id
    mat = results[-1].get_material(str(mat_id))
    sources[uid] = mat.decay_photon_source

with open('sources.pkl', 'wb') as fh:
    dill.dump(sources, fh)
