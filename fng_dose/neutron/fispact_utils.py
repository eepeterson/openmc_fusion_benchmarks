import openmc
import numpy as np

def fispact_photon_energy(timestep, density):
    spectrum = timestep.gamma_spectrum
    energies_ev = np.asarray(spectrum.boundaries) * 1e6
    cm3 = timestep.total_mass * 1e3 / density
    gammas_per_sec = np.asarray(spectrum.volumetric_rates) * cm3
    gammas_per_ev_sec = gammas_per_sec / np.diff(energies_ev)
    gammas_per_ev_sec = np.concatenate((gammas_per_ev_sec, [0.0]))
    return openmc.stats.Tabular(energies_ev, gammas_per_ev_sec, 'histogram')


def fispact_material(timestep, density):
    mat = openmc.Material()
    mat.volume = timestep.total_mass * 1e3 / density
    for nuc in timestep.nuclides:
        if 'm' in nuc.state:
            name = f'{nuc.element}{nuc.isotope}_m1'
        else:
            name = f'{nuc.element}{nuc.isotope}'
        mat.add_nuclide(name, nuc.atoms / mat.volume * 1.0e-24)
    mat.set_density('g/cm3', density)
    return mat
