# R2S Calculation

## Preprocessing steps

- Run `r2s.py --run-volume-calc` to calculate volume for each activation cell
  (written to cell_volumes.json)
- Run `reduce_chain.py` to create a reduced depletion chain (in the future,
  this shouldn't be necessary but right now the `reduce_chain` argument doesn't
  work when C0 is present)

## Full R2S calculation

Run `r2s.py` to run the neutron transport, activation, and photon transport
calculation. An intermediate `activation/sources.pkl` file is create with the
generate decay photon source. This script has an `--operator` command-line
argument that can be set to "coupled" or "independent", controlling which
transport operator is used in the calculation.
