# R2S Calculation

## Preprocessing steps

1. Run `generate_settings.py` to produce neutron source in settings.xml
2. Run `volumes.py` to calculate volume for each activation cell
3. Run `reduce_chain.py` to create a reduced depletion chain (in the future,
   this shouldn't be necessary but right now the `reduce_chain` argument doesn't
   work when C0 is present)

## Full R2S calculation

Run `r2s.py <cooling_time>` to run the neutron transport, activation, and photon
transport calculation. An intermediate `activation/sources.pkl` file is create
with the generate decay photon source. This script has an `--operator`
command-line argument that can be set to "coupled" or "independent", controlling
which transport operator is used in the calculation.
