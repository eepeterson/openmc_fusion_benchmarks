[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)
[![CI testing](https://github.com/SteSeg/openmc_fusion_benchmarks/actions/workflows/ci.yml/badge.svg?branch=add_ci)](https://github.com/SteSeg/openmc_fusion_benchmarks/workflows/ci.yml)
[![Code Coverage](https://coveralls.io/repos/github/SteSeg/openmc_fusion_benchmarks/badge.svg?branch=add_ci)](https://coveralls.io/github/SteSeg/openmc_fusion_benchmarks?branch=add_ci)
[![Jupyter Book Badge](https://jupyterbook.org/badge.svg)](https://eepeterson.github.io/openmc_fusion_benchmarks)



# OpenMC Fusion Benchmarks

**Next-generation CAD-based V&V repository** for fusion neutronics with schema-validated, code-agnostic benchmark definitions.

## Key Features

- ✅ **Schema-validated specifications** - Every benchmark validated against [`benchmark_schema.yaml`](src/openmc_fusion_benchmarks/benchmarks/benchmark_schema.yaml)
- 🔧 **CAD-based geometry** - Native CAD files automatically meshed to DAGMC
- 🔄 **Code-agnostic** - Same `specifications.yaml` → reproducible across any transport code
- 📊 **Built-in results database** - Experimental and computational results packaged with the repository

## Installation

```bash
git clone --recurse-submodules https://github.com/eepeterson/openmc_fusion_benchmarks.git
cd openmc_fusion_benchmarks
pip install .
```

Requirements: Python ≥3.7, OpenMC ≥0.14.0

## Quick Start

**Run a benchmark:**
```python
import openmc_fusion_benchmarks as ofb

benchmark = ofb.OpenmcBenchmark(name="oktavian_al")
benchmark.run()
```

## Documentation

[Read the Docs](https://openmc-fusion-benchmarks.readthedocs.io/en/latest/index.html)

## License

MIT License
