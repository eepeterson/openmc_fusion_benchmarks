from .sandy_wrapper import perturb_to_hdf5
from .data_conventions import get_nuclide_gnds
from .modify_xs_xml import rewrite_xs_xml, perturb_xs_xml
import os
import openmc
import numpy as np


def tmc_engine(model: openmc.Model, nsamples: int, lib_name: str, nuclide,
               reaction: int, perturb_xs: bool = False):
    # convert nuclide to gnds name
    nuclide = get_nuclide_gnds(nuclide)
    xs_file = f'cross_sections_mod.xml'
    path_to_file = f'results_{nuclide}.h5'

    if perturb_xs:
        perturb_to_hdf5(nsamples, lib_name, nuclide, reaction, nprocesses=1,
                        error=.001)

    for n in np.arange(nsamples):

        rewrite_xs_xml(new_xs_file=xs_file)

        directory = f"{nuclide}_{lib_name}"
        xs_h5_file = f"{directory}/{nuclide}_{n}_{lib_name}.h5"
        perturb_xs_xml(xs_file, xs_h5_file, nuclide)
        openmc.config['cross_sections'] = xs_file

        # run simulation
        model.run()

        # postprocess result
        sp = openmc.StatePoint('statepoint.100.h5')
        # open tally and push to hdf5
        for t in sp.tallies:
            tally = sp.get_tally(id=t)
            tally_df = tally.get_pandas_dataframe()
            tally_df.to_hdf(path_to_file, tally.name+f'_{n}', mode='a',
                            format='table', data_columns=True, index=False)

        os.remove('summary.h5')
        os.remove('statepoint.100.h5')
