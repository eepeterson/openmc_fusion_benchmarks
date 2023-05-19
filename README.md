# openmc_sinbad_benchmarks

## Installation

Before starting running the benchmarks install the `benchmark_tools` package. From the cloned repo:
```
cd benchmark_tools
pip install .
```

## Folder organization

- The `benchmark_tools` folder is a package that contains some routines that are imported in most of the models. For instance the `fng_source` is a piece of code that models the Frascati Neutron Generator (FNG) source according to `OpenMC` sintax. While, the `From_IRDFF` routine allows to include the dosimetry cross sections often used for modeling activation foils.
- The `IRDFF2_xs` folder contains the dosimetry cross sections required for modeling the activation foils in the models
- The `models` folder contains subfolders with all the benchmarks modeled so far with `OpenMC`
