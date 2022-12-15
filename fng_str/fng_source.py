# %%
# building the FNG neutron source
# angular and energy distribution

import pandas as pd
import numpy as np
import openmc


def fng_source(center=(0, 0, 0), reference_uvw=(1, 0, 0)):
    '''method for building the Frascati Neutron Generator source in OpenMC
    
    Parameters
    ----------
    center : coordinate position of the source (it is a point source)

    reference_uvw : direction for the polar angle (tuple or list of versors)
    it is the same for the openmc.PolarAzimuthal class
    more specifically, polar angle = 0 is the direction of the D accelerator
    towards the Ti-T target    
    '''

    # read tabulated data
    source_energy = pd.read_excel(r'fng_source.xlsx', sheet_name='Energy')
    source_yield = pd.read_excel(r'fng_source.xlsx', sheet_name='Yield')
    energies = [source_energy[column][1:]*1e6 for column in source_energy]
    yields = [source_yield[column][1:] for column in source_yield]

    source_strength = pd.read_excel(r'fng_source.xlsx', sheet_name='rel_strength')
    tot_strength = sum(source_strength['Relative Intensity'])

    # angles for the angular distribution
    theta = np.linspace(0, 180, 19)
    polar_bounds = np.cos(np.radians(np.r_[0, theta[:-1] + 5, 180]))

    # azimuthal angle
    phi = openmc.stats.Uniform(a=0, b=2*np.pi)

    # generating sources
    all_sources = []
    for ab in range(len(theta)):

        mu = openmc.stats.Uniform(a=polar_bounds[ab], b=polar_bounds[ab+1])
        space = openmc.stats.Point(center)

        # # cylindrical geometry
        # r_cyl = openmc.stats.Uniform(a=0, b=1.2)
        # phi_cyl = openmc.stats.Uniform(a=0, b=2*np.pi)
        # z_cyl = openmc.stats.Uniform(a=0, b=1.0)
        # origin_cyl = (0, -0.5, 0)
        # space = openmc.stats.CylindricalIndependent(r=r_cyl, phi=phi_cyl, z=z_cyl, origin=origin_cyl)
        # #

        angle = openmc.stats.PolarAzimuthal(mu=mu, phi=phi, reference_uvw=reference_uvw)
        energy = openmc.stats.Tabular(energies[ab], yields[ab])
        strength = source_strength['Relative Intensity'][ab] / tot_strength

        source = openmc.Source(space=space, angle=angle, energy=energy, strength=strength, particle='neutron')

        all_sources.append(source)
    
    return all_sources


# %%



