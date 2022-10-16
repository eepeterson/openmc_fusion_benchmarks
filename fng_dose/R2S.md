# R2S Calculation

## Neutron transport + activation

From the `neutron/` directory:

1. Run `generate_settings.py` to produce neutron source in settings.xml
2. Run `add_volumes.py` to add volume information to each material
3. Run `reduce_chain.py` to create a reduced depletion chain (in the future,
   this shouldn't be necessary but right now the `reduce_chain` argument doesn't
   work when C0 is present)
4. Run `activate.py` to run the activation calculation and create sources.pkl

## Photon transport

From the `photon/` directory:

1. Run `run_photon_dose.py <time>` to generate photon sources based on
   the sources.pkl file at the specified cooling time.
