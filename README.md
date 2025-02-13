[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![CI testing](https://github.com/SteSeg/openmc_fusion_benchmarks/workflows/ci.yml/badge.svg?branch=main)](https://github.com/SteSeg/openmc_fusion_benchmarks/workflows/ci.yml)
[![cov](https://SteSeg.github.io/openmc_fusion_benchmarks/badges/coverage.svg)](https://github.com/SteSeg/openmc_fusion_benchmarks/actions)

# openmc_fusion_benchmarks
OpenMC Fusion Benchmarks is a platform for V&V of fusion neutronics. It focuses on fusion-relevant integral benchmarks. It relies on an automated workflow for model simulation, data postprocessing, visualization and analysis. It embeds a database of experimental and numerical results for quick comparisons to which users can contribute to. The database contribution workflow is fully automated.

## Installation
```
pip install -e .
```

## Folder organization
The `src/` folder provides all the classes and functions of the package. It also has additional data in the `data/` subfolder (e.g. IRDFF-2 cross sections) and some neutron source characteristics, importable as `openmc.Source` object, in the `neutron_sources/` subfolder (so far there is only the Frascati Neutron Generator - FNG energy and angular characteristics).

The `notebooks/` folder contains some useful tutorial notebooks.

The `models/` folder has all the models and benchmarks covered by this repo (e.g. the `fng-str` and `fns-duct` benchmarks). Each bechmark provides an `openmc_model.py` file that has the Python API input file for the relative openmc model. It can be run via command line and often has some different options for running it. There should be also a `postprocessing.ipynb` jupyter notebook for quick result visualization and confrontation against other experimental or other codes' results. The `run_and_store.py` script is an automated way to run the openmc model, postprocess the results and store them in a hdf5 file in the `results_database/` in just one click. It sometimes runs the model more than once in order match all the results provided by the experiments (e.g. the `fng-str` experiment has three different configurations and thus requires three different simulations).

## Interactive database
Each model has a `results_database/` subfolder that contains the experimental results as well as previous simulation results of the benchmark in form of hdf files. Tutorial notebooks in the `notebooks/` folder show how to navigate those files. All the hdf files belonging to the same model have same structure (same tallies with the same name, postprocessed in the same way and extractable in the same way). Each `results_database/` subfolder should have at least the `experiment.h5` file. Results from openmc and other codes are more than welcome, especially if built consistently with the files already present. The result hdf5 files from codes should have the name of `codename-codeversion_xsnameandversion.h5`.

## Database contribution
The contribution workflow is automated for the openmc code. Users can run any openmc model and then the postprocessing notebook for a direct validation of the openmc version they use but they can also contribute to the repository. To contribute it is sufficient to:

- Run the `run_and_store.py` script of a benchmark of choice according to the `run_and_store.ipynb` tutorial notebook (this would generate a result hdf5 file in the respective `results_database/` folder)
- Review the result hdf5 and open a pull request to propose the new file to the database
