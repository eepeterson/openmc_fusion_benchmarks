# Developer Installation

First clone the repository

```bash
git clone git@github.com:eepeterson/openmc_fusion_benchmarks.git
```

The change directory into the repository

```bash
cd openmc_fusion_benchmarks
```

Then install the package including the dependencies needed for testing and building the documentation

```bash
python -m pip install -e .[docs,tests]
```