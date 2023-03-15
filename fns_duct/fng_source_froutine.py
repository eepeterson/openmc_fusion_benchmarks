# %%
# building the FNG neutron source
# angular and energy distribution

import numpy as np
import openmc


def fng_source(center=(0, 0, 0), reference_uvw=(0, 0, 1)):
    '''method for building the Frascati Neutron Generator source in OpenMC
    with data tabulated from the fortran->C++ routine
    
    Parameters
    ----------
    center : coordinate position of the source (it is a point source)

    reference_uvw : direction for the polar angle (tuple or list of versors)
    it is the same for the openmc.PolarAzimuthal class
    more specifically, polar angle = 0 is the direction of the D accelerator
    towards the Ti-T target    
    '''

    # read tabulated data
    fng_source_fr = np.loadtxt("fng_from_routine.csv", delimiter=",")

    # angular bins in [0, pi)
    pbins = np.cos(np.linspace(0, np.pi, 37))

    # energy and flux values from tables
    evalues = fng_source_fr[0]
    fvalues = fng_source_fr[1:]

    # yield values for strengths
    yields = np.sum(fvalues, axis=-1) * np.diff(pbins)
    yields /= np.sum(yields)

    # azimuthal values
    phi = openmc.stats.Uniform(a=0, b=2*np.pi)

    all_sources = []
    for i,angle in enumerate(pbins[:-1]):

        mu = openmc.stats.Uniform(a=pbins[i+1], b=pbins[i])

        space = openmc.stats.Point(center)
        angle = openmc.stats.PolarAzimuthal(mu=mu, phi=phi, reference_uvw=reference_uvw)
        energy = openmc.stats.Tabular(evalues, fvalues[i], interpolation='linear-linear')
        strength = yields[i]

        my_source = openmc.Source(space=space, angle=angle, energy=energy, strength=strength, particle='neutron')

        all_sources.append(my_source)

    
    return all_sources


# %%



