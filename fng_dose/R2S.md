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

1. Run `generate_sources_tallies.py` to generate photon sources based on the
   sources.pkl file
2. Run `apply_source_cutoff.py E` to restrict the source to only photons above a
   certain energy.
2. Run `openmc`
3. Run `evaluate_dose.py` to calculate the dose in Sv/h

Right now step 2 is necessary because with TTB on, the minimum photon energy is
slightly less than 1100.0 eV. Restricting to 1100.0 eV allows OpenMC to run.
