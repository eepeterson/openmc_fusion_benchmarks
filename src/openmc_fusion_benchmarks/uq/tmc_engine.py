from pathlib import Path

import numpy as np
import openmc
import xarray as xr

from .uq_utils import perturb_xs_xml, get_nuclide_gnds, perturb_to_hdf5
from ..utils import _openmc_to_ofb, _save_result


def tmc_engine(model: openmc.Model, realizations: int, lib_name: str, nuclide,
               reaction: int = None, perturb_xs: bool = True,
               _is_benchmark: bool = False, _mesh=None, _spec_tallies=None,
               *args, **kwargs):
    """Runs a TMC simulation on a given OpenMC model object. 
    With perturb_xs=True it is possible to perturb the cross sections of a 
    specific nuclide and reaction from a given nuclear data library 
    automatically before the starting of the actual TMC simulation.
    The results of the TMC simulation are stored in a .h5 file as OpenMC
    tallies in Xarray DataArray format.

    Parameters
    ----------
    model : openmc.Model
        OpenMC model object to run TMC simulations on
    realizations : int
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

    # runs sandy and generates perturbed xs only if perturb_xs is True
    if perturb_xs:
        perturb_to_hdf5(realizations, lib_name, nuclide, reaction, nprocesses=1,
                        error=.001)

    for n in np.arange(realizations):

        directory = f"{nuclide}_{lib_name}"
        xs_h5_file = f"{directory}/{nuclide}_{n}_{lib_name}.h5"
        perturb_xs_xml(xs_file, xs_h5_file, nuclide)
        openmc.config['cross_sections'] = xs_file

        # run simulation
        statepoint = model.run(*args, **kwargs)

        # Postprocess results in case of standalone TMC run
        sp = openmc.StatePoint(statepoint)

        # realization_label = f'{nuclide}_{n}_{lib_name}'
        realization_label = f'{nuclide}_{n}'
        if _is_benchmark:
            # Postprocess results in case of OFB benchmark model
            _openmc_to_ofb(statepoint=sp,
                           realization_label=realization_label,
                           mesh=_mesh,
                           spec_tallies=_spec_tallies,
                           *args, **kwargs)
        else:
            # open tally and push to netcdf ofb style
            for t in sp.tallies:
                tally = sp.get_tally(id=t)
                df = tally.get_pandas_dataframe()

                df = df.drop(columns=['surface', 'cell', 'particle', 'nuclide',
                              'score', 'energyfunction'], errors='ignore')

                # Convert to xarray and add dimensions
                da = xr.DataArray(
                    df.values[np.newaxis, :, :],  # shape: (1, r, c)
                    dims=["realization", "row", "column"],
                    coords={
                        "realization": [realization_label],
                        "column": df.columns,
                        "row": np.arange(df.shape[0]),
                    },
                    name=tally.name
                )

                _save_result(new_result=da, filename="tmc_results.h5",
                             group=tally.name, realization_label=f'{nuclide}_{n}_{lib_name}')

        Path('summary.h5').unlink(missing_ok=True)
        Path(statepoint).unlink(missing_ok=True)
