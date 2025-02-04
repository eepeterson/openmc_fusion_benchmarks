from .sandy_wrapper import perturb_to_hdf5
from .data_conventions import get_nuclide_gnds
from .modify_xs_xml import rewrite_xs_xml, perturb_xs_xml
import os
import openmc
import numpy as np


def tmc_engine(model: openmc.Model, nsamples: int, lib_name: str, nuclide,
               reaction: int = None, perturb_xs: bool = True):
    """Runs a TMC simulation on a given OpenMC model object. 
    With perturb_xs=True it is possible to perturb the cross sections of a 
    specific nuclide and reaction from a given nuclear data library 
    automatically before the starting of the actual TMC simulation.
    The results of the TMC simulation are stored in a .h5 file as OpenMC
    tallies in PandasDataFrame format.

    Parameters
    ----------
    model : openmc.Model
        OpenMC model object to run TMC simulations on
    nsamples : int
        Number of samples to run in the TMC simulation 
        (i.e. number of times the xs is perturbed)
    lib_name : str
        Name of the nuclear data library to perturb the xs from 
        (e.g. 'ENDF/B-VIII.0')
    nuclide : str or int
        Identifier of the nuclide for which the cross section is perturbed
        (GNDS, ZAID or ZAM)
    reaction : int, optional
        MT value for the specific reaction to perturb, by default None
    perturb_xs : bool, optional
        Flag for the automatic generation of perturbed .h5 xs files right
        before running the TMC simulation. Set to False if the use has already
        a set of perturbed xs in .h5 format
        to point to, by default True
    """

    # convert nuclide to gnds name
    nuclide = get_nuclide_gnds(nuclide)
    xs_file = f'cross_sections_mod.xml'
    path_to_file = f'tmc_results_{nuclide}.h5'

    # runs sandy and generates perturbed xs only if perturb_xs is True
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
        sp_name = f'statepoint.{model.settings.batches}.h5'
        sp = openmc.StatePoint(sp_name)
        # open tally and push to hdf5
        for t in sp.tallies:
            tally = sp.get_tally(id=t)
            tally_df = tally.get_pandas_dataframe()
            tally_df.to_hdf(path_to_file, tally.name+f'_{n}', mode='a',
                            format='table', data_columns=True, index=False)

        os.remove('summary.h5')
        os.remove(sp_name)
