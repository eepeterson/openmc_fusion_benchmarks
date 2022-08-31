from math import pi
from uncertainties import ufloat
import openmc

# Get volume of cell 651 (radius of 1.9 cm)
r = 1.9
cm3 = 4/3 * pi * r**3
g_per_cm3 = 0.0011050927231898923
kg = g_per_cm3*cm3 * 1e-3

with openmc.StatePoint('statepoint.40.h5') as sp:
    tallies = list(sp.tallies.values())
    flux_tally, heating_tally = tallies

# Calculate dose using flux-to-dose conversion factor
mean = flux_tally.mean.ravel()[0]
stdev = flux_tally.std_dev.ravel()[0]
mrem_cm3_per_sec = ufloat(mean, stdev)
mrem_per_sec = mrem_cm3_per_sec / cm3
Sv_per_mrem = 1e-5
Sv_per_h = mrem_per_sec * 3600. * Sv_per_mrem
print(f'Dose rate (flux) = {Sv_per_h} Sv/h')

# Calculate dose using energy deposition directly
mean = heating_tally.mean.ravel()[0]
stdev = heating_tally.std_dev.ravel()[0]
eV_per_sec = ufloat(mean, stdev)
J_per_eV = 1.602176634e-19
J_per_kg_per_sec = eV_per_sec / kg * J_per_eV
Sv_per_h = J_per_kg_per_sec * 3600.
print(f'Dose rate (Edep) = {Sv_per_h} Sv/h')
