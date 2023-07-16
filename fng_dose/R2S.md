# R2S Calculation

## Preprocessing steps

1. Run `generate_settings.py` to produce neutron source in settings.xml
2. Run `volumes.py` to calculate volume for each activation cell
3. Run `reduce_chain.py` to create a reduced depletion chain (in the future,
   this shouldn't be necessary but right now the `reduce_chain` argument doesn't
   work when C0 is present)

## Neutron transport + activation

Run `activate.py` to run the activation calculation and create
`activation/sources.pkl`. This script has an `--operator` command-line argument
that can be set to "coupled" or "independent", controlling which transport
operator is used in the calculation.

## Photon transport

Run `run_photon_dose.py <time>` to generate photon sources based on the
`activation/sources.pkl` file at the specified cooling time.
