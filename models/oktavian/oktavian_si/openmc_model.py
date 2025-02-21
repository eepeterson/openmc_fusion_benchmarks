#!/usr/bin/env python3
import argparse
import numpy as np

import openmc


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batches", type=int, default=100,
                        help='Number of batches to simulate (int)')
    parser.add_argument("-p", "--particles", type=int, default=int(1e7),
                        help='Number of particles per batch (int)')
    parser.add_argument("-s", "--threads", type=int,
                        help='Number of threads to use in the simulation (int)')

    args = parser.parse_args()

    return args

def main():
    """Analysis of Osaka Sphere Benchmark Experiment (OKTAVIAN) -Tungsten experiment"""
    
    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    ############################################################################
    # Define Materials

    # 
    mat_1 = openmc.Material(material_id=1)
    mat_1.set_density('g/cm3', 1.29581)
    mat_1.add_nuclide('Si28', 0.9223, 'ao')
    mat_1.add_nuclide('Si29', 0.04683, 'ao')
    mat_1.add_nuclide('Si30', 0.03087, 'ao')
    # 
    mat_2 = openmc.Material(material_id=2)
    mat_2.set_density('g/cm3', 7.824)
    mat_2.add_nuclide('Cr50', 0.00803825, 'wo')
    mat_2.add_nuclide('Cr52', 0.15501, 'wo')
    mat_2.add_nuclide('Cr53', 0.0175768, 'wo')
    mat_2.add_nuclide('Cr54', 0.00437525, 'wo')
    mat_2.add_nuclide('Fe54', 0.0411488, 'wo')
    mat_2.add_nuclide('Fe56', 0.645948, 'wo')
    mat_2.add_nuclide('Fe57', 0.0149178, 'wo')
    mat_2.add_nuclide('Fe58', 0.00198528, 'wo')
    mat_2.add_nuclide('Ni58', 0.0755655, 'wo')
    mat_2.add_nuclide('Ni60', 0.0291075, 'wo')
    mat_2.add_nuclide('Ni61', 0.0012654, 'wo')
    mat_2.add_nuclide('Ni62', 0.00403374, 'wo')
    mat_2.add_nuclide('Ni64', 0.00102786, 'wo')
    
    # create materials instance
    materials = openmc.Materials([mat_1, mat_2])

    # %%

    ############################################################################
    # Define Geometry

    # surfaces
    surf_3 = openmc.Sphere(surface_id=3, x0=0.0, y0=0.0, z0=0.0, r=10.0)
    surf_8 = openmc.XPlane(surface_id=8, x0=8.32)
    surf_1 = openmc.XCylinder(surface_id=1, y0=0.0, z0=0.0, r=5.55)
    surf_6 = openmc.Sphere(surface_id=6, x0=0.0, y0=0.0, z0=0.0, r=30.5)
    surf_4 = openmc.Sphere(surface_id=4, x0=0.0, y0=0.0, z0=0.0, r=10.2)
    surf_2 = openmc.XCylinder(surface_id=2, y0=0.0, z0=0.0, r=5.75)
    surf_5 = openmc.Sphere(surface_id=5, x0=0.0, y0=0.0, z0=0.0, r=30.0)
    surf_7 = openmc.Sphere(surface_id=7, x0=0.0, y0=0.0, z0=0.0, r=100.0, boundary_type='vacuum')

    # regions
    region_1 = ((-surf_3 & -surf_8) | (+surf_8 & -surf_1 & -surf_6))
    region_2 = ((+surf_3 & -surf_4 & -surf_8) | (+surf_8 & +surf_1 & -surf_2 & -surf_6))
    region_3 = ((+surf_4 & -surf_5 & -surf_8) | (+surf_8 & +surf_2 & -surf_5))
    region_4 = ((+surf_5 & -surf_6 & -surf_8) | (+surf_8 & +surf_2 & +surf_5 & -surf_6))
    region_5 = (+surf_6 & -surf_7)

    # cells
    cell_1 = openmc.Cell(cell_id=1, region=region_1, fill=None)
    cell_2 = openmc.Cell(cell_id=2, region=region_2, fill=mat_2)
    cell_3 = openmc.Cell(cell_id=3, region=region_3, fill=mat_1)
    cell_4 = openmc.Cell(cell_id=4, region=region_4, fill=mat_2)
    cell_5 = openmc.Cell(cell_id=5, region=region_5, fill=None)

    # create root universe
    universe = openmc.Universe(cells=[cell_1, cell_2, cell_3, cell_4, cell_5])

    # create geometry instance
    model.geometry = openmc.Geometry(universe)

    ############################################################################
    # Define Settings

    # source definition
    # oktavian source
    # Create a point source
    energies = np.array([
        0.14489, 0.15696, 0.17003, 0.18419, 0.19953, 0.21615, 0.23415, 0.25365, 
        0.27478, 0.29766, 0.32245, 0.34931, 0.37840, 0.40992, 0.44406, 0.48104, 
        0.52111, 0.56451, 0.61153, 0.66246, 0.71763, 0.77740, 0.84215, 0.91229, 
        0.98827, 1.07060, 1.15980, 1.25630, 1.36100, 1.47430, 1.59710, 1.73010, 
        1.87420, 2.03030, 2.19940, 2.38260, 2.58110, 2.79600, 3.02890, 3.28120, 
        3.55450, 3.85050, 4.17120, 4.51860, 4.89490, 5.30260, 5.74430, 6.22270, 
        6.74100, 7.30240, 7.91060, 8.56940, 9.28310, 10.05600, 10.89400, 11.80100, 
        12.78400, 13.84900, 15.00200, 16.25200, 17.60500, 19.07200, 20.66000
    ])*1e6

    weights = np.array([
        5.526160e-05, 4.505010e-05, 7.538530e-05, 0.0001525780, 0.0003242640, 
        0.0001333580, 0.0002099350, 0.0002942890, 0.0003575800, 0.0003455120, 
        0.0002985720, 0.0004691770, 0.0004705320, 0.0005297930, 0.0005564800, 
        0.0000665921, 0.0006847820, 0.0007484870, 0.0008127610, 0.0009074510, 
        0.0009113700, 0.0010903600, 0.0011200000, 0.0012695700, 0.0011870000, 
        0.0014949300, 0.0015342100, 0.0015147600, 0.0015700800, 0.0015923600, 
        0.0015356400, 0.0015008500, 0.0016821100, 0.0016157800, 0.0016312000, 
        0.0016782500, 0.0030940400, 0.0028755500, 0.0014989800, 0.0014188900, 
        0.0014437900, 0.0013049700, 0.0013850900, 0.0012661900, 0.0011862100, 
        0.0010310200, 0.0011861000, 0.0010952600, 0.0010458500, 0.0012784900, 
        0.0013835800, 0.0015852400, 0.0017453500, 0.0024926900, 0.0034166500, 
        0.0059493100, 0.0158133000, 0.6004580000, 0.3119180000, 0.0063896000, 
        8.341600e-05, 4.769740e-05, 0.0
    ])

    source = openmc.IndependentSource()
    source.space = openmc.stats.Point((0, 0, 0))
    source.angle = openmc.stats.Isotropic()
    source.particle = 'neutron' 
    source.energy = openmc.stats.Tabular(energies, weights, interpolation='histogram')
    
    # Indicate how many particles to run
    settings = openmc.Settings(run_mode='fixed source')
    settings.photon_transport = True
    settings.source = source
    settings.batches = args.batches
    settings.particles = args.particles
    settings.output = {'tallies': False}

    model.settings = settings

    ############################################################################
    # Specify Tallies
    # setup the filters for the cell tally
    vessel_surface_filter = openmc.SurfaceFilter(surf_6) # outer radius of outer steel vessel
    neutron_particle_filter = openmc.ParticleFilter(['neutron'])
    photon_particle_filter = openmc.ParticleFilter(['photon'])

    neutron_energies = np.array([
        9.14600E-02, 9.61500E-02, 1.01080E-01, 1.06260E-01, 1.11710E-01, 1.17440E-01,
        1.23460E-01, 1.29790E-01, 1.36440E-01, 1.43440E-01, 1.50790E-01, 1.58530E-01,
        1.66650E-01, 1.75200E-01, 1.84180E-01, 1.93620E-01, 2.03550E-01, 2.13990E-01,
        2.24960E-01, 2.36490E-01, 2.48620E-01, 2.61360E-01, 2.74760E-01, 2.88850E-01,
        3.03660E-01, 3.19230E-01, 3.35600E-01, 3.52800E-01, 3.70890E-01, 3.89910E-01,
        4.09900E-01, 4.30920E-01, 4.53010E-01, 4.76240E-01, 5.00650E-01, 5.26320E-01,
        5.53310E-01, 5.81680E-01, 6.11500E-01, 6.42850E-01, 6.75810E-01, 7.10460E-01,
        7.46890E-01, 7.85180E-01, 8.25440E-01, 8.67760E-01, 9.12250E-01, 9.59020E-01,
        1.00820E+00, 1.05990E+00, 1.11420E+00, 1.17140E+00, 1.23140E+00, 1.29450E+00,
        1.36090E+00, 1.43070E+00, 1.50400E+00, 1.58120E+00, 1.66220E+00, 1.74750E+00,
        1.83700E+00, 1.93120E+00, 2.03020E+00, 2.13430E+00, 2.24380E+00, 2.35880E+00,
        2.47980E+00, 2.60690E+00, 2.74050E+00, 2.88110E+00, 3.02880E+00, 3.18410E+00,
        3.34730E+00, 3.51890E+00, 3.69930E+00, 3.88900E+00, 4.08840E+00, 4.29900E+00,
        4.51840E+00, 4.75010E+00, 4.99360E+00, 5.24960E+00, 5.51880E+00, 5.80170E+00,
        6.09920E+00, 6.41190E+00, 6.74060E+00, 7.08620E+00, 7.44960E+00, 7.83150E+00,
        8.23300E+00, 8.65520E+00, 9.09890E+00, 9.56540E+00, 1.00560E+01, 1.05710E+01,
        1.11130E+01, 1.16830E+01, 1.22820E+01, 1.29120E+01, 1.35740E+01, 1.42700E+01,
        1.50020E+01, 1.57710E+01, 1.65790E+01, 1.74290E+01, 1.83230E+01, 1.92620E+01,
        2.02500E+01
    ])*1e6

    photon_energies = np.array([
        3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 7.00000E-01, 8.00000E-01,
        9.00000E-01, 1.00000E+00, 1.10000E+00, 1.20000E+00, 1.30000E+00, 1.40000E+00,
        1.50000E+00, 1.60000E+00, 1.70000E+00, 1.80000E+00, 1.90000E+00, 2.00000E+00,
        2.10000E+00, 2.20000E+00, 2.30000E+00, 2.40000E+00, 2.50000E+00, 2.60000E+00,
        2.70000E+00, 2.80000E+00, 2.90000E+00, 3.00000E+00, 3.10000E+00, 3.20000E+00,
        3.30000E+00, 3.40000E+00, 3.50000E+00, 3.60000E+00, 3.70000E+00, 3.80000E+00,
        3.90000E+00, 4.00000E+00, 4.50000E+00, 5.00000E+00, 5.50000E+00, 6.00000E+00,
        6.50000E+00, 7.00000E+00, 7.50000E+00, 8.00000E+00, 8.50000E+00, 9.00000E+00,
        9.50000E+00, 1.00000E+01, 1.05000E+01, 1.10000E+01
    ])*1e6

    energy_neutron_filter = openmc.EnergyFilter(neutron_energies)
    energy_photon_filter = openmc.EnergyFilter(photon_energies)

    #creates an empty tally object
    model.tallies = openmc.Tallies()

    # create the tally
    neutron_surface_spectra_tally = openmc.Tally(name='nspectrum')
    neutron_surface_spectra_tally.scores = ['current']
    neutron_surface_spectra_tally.filters = [vessel_surface_filter, neutron_particle_filter, energy_neutron_filter]
    model.tallies.append(neutron_surface_spectra_tally)

    photon_surface_spectra_tally = openmc.Tally(name='gspectrum')
    photon_surface_spectra_tally.scores = ['current']
    photon_surface_spectra_tally.filters = [vessel_surface_filter, photon_particle_filter, energy_photon_filter]
    model.tallies.append(photon_surface_spectra_tally)

    # define the folder names for storing the statepoints
    cwd = 'results'
    model.settings = settings
    
    return model.run(cwd=cwd, threads=args.threads)


if __name__ == "__main__":
    main()