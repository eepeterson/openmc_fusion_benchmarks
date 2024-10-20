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
    mat_1.set_density('g/cm3', 3.72)
    mat_1.add_nuclide('Cr50', 0.0433683, 'wo')
    mat_1.add_nuclide('Cr52', 0.836315, 'wo')
    mat_1.add_nuclide('Cr53', 0.0948314, 'wo')
    mat_1.add_nuclide('Cr54', 0.0236055, 'wo')
    mat_1.add_nuclide('Fe54', 9.352e-05, 'wo')
    mat_1.add_nuclide('Fe56', 0.00146806, 'wo')
    mat_1.add_nuclide('Fe57', 3.3904e-05, 'wo')
    mat_1.add_nuclide('Fe58', 4.512e-06, 'wo')
    mat_1.add_nuclide('C12', 0.0002074814374679942, 'wo')
    mat_1.add_nuclide('C13', 2.5185625320058296e-06, 'wo')
    mat_1.add_nuclide('Si28', 6.4561e-05, 'wo')
    mat_1.add_nuclide('Si29', 3.2781e-06, 'wo')
    mat_1.add_nuclide('Si30', 2.1609e-06, 'wo')
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
    surf_6 = openmc.Sphere(surface_id=6, x0=0.0, y0=0.0, z0=0.0, r=19.95)
    surf_4 = openmc.Sphere(surface_id=4, x0=0.0, y0=0.0, z0=0.0, r=10.2)
    surf_2 = openmc.XCylinder(surface_id=2, y0=0.0, z0=0.0, r=5.75)
    surf_5 = openmc.Sphere(surface_id=5, x0=0.0, y0=0.0, z0=0.0, r=19.75)
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
        0.10000, 0.10417, 0.11172, 0.11982, 0.12851, 0.13783, 0.14782, 0.15854, 
        0.17003, 0.18236, 0.19559, 0.20977, 0.22498, 0.24129, 0.25878, 0.27755, 
        0.29767, 0.31926, 0.34241, 0.36723, 0.39386, 0.42242, 0.45305, 0.48590, 
        0.52113, 0.55891, 0.59944, 0.64290, 0.68952, 0.73952, 0.79314, 0.85064, 
        0.91232, 0.97847, 1.04940, 1.12550, 1.20710, 1.29470, 1.38850, 1.48920, 
        1.59720, 1.71300, 1.83720, 1.97040, 2.11330, 2.26650, 2.43080, 2.60710, 
        2.79610, 2.99890, 3.21630, 3.44950, 3.69970, 3.96790, 4.25560, 4.56420, 
        4.89510, 5.25010, 5.63070, 6.03900, 6.47690, 6.94650, 7.45020, 7.99040, 
        8.56970, 9.19110, 9.85760, 10.57200, 11.33900, 12.16100, 13.04300, 13.98900, 
        15.00300, 16.09100, 17.25700, 18.50900, 19.85100, 21.29000
    ])*1e6

    weights = np.array([
        0.001404640, 0.000000000, 0.002002340, 0.000906709, 0.000363535,
        0.000317011, 0.000279367, 1.82539e-05, 0.000209078, 0.000343060, 
        0.000379611, 0.000304170, 0.000303979, 0.000350817, 0.000393017, 
        0.000400549, 0.000418799, 0.000434474, 0.000515476, 0.000642848, 
        0.000652009, 0.000634607, 0.000738784, 0.000757084, 0.000839447, 
        0.000938240, 0.001036530, 0.001056820, 0.001186460, 0.001220910, 
        0.001273590, 0.001327930, 0.001426100, 0.001401420, 0.001528710, 
        0.001672700, 0.001670330, 0.001630910, 0.001688410, 0.001688710, 
        0.001723820, 0.001710920, 0.001709670, 0.001724100, 0.001751820, 
        0.001737340, 0.001774330, 0.001904000, 0.001770590, 0.001624660, 
        0.001641450, 0.001579430, 0.001486700, 0.001459820, 0.001433400, 
        0.001445320, 0.001418450, 0.001289990, 0.001258470, 0.001247670, 
        0.001287230, 0.001238070, 0.001291640, 0.001505670, 0.002030290, 
        0.002325010, 0.003062420, 0.004728200, 0.007763550, 0.016588500, 
        0.082441200, 0.711895000, 0.105764000, 0.001808060, 0.000439258, 
        6.24285e-05, 3.95600e-05, 0.0
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
        1.00080E-01, 1.03120E-01, 1.06260E-01, 1.09500E-01, 1.12830E-01, 1.16270E-01,
        1.19810E-01, 1.23460E-01, 1.27220E-01, 1.31090E-01, 1.35090E-01, 1.39200E-01,
        1.43440E-01, 1.47810E-01, 1.52310E-01, 1.56950E-01, 1.61730E-01, 1.66650E-01, 
        1.71730E-01, 1.76960E-01, 1.82350E-01, 1.87900E-01, 1.93620E-01, 1.99520E-01, 
        2.05600E-01, 2.11860E-01, 2.18310E-01, 2.24960E-01, 2.31810E-01, 2.38870E-01, 
        2.46140E-01, 2.53640E-01, 2.61370E-01, 2.69320E-01, 2.77530E-01, 2.85980E-01,
        2.94690E-01, 3.03660E-01, 3.12910E-01, 3.22440E-01, 3.32260E-01, 3.42380E-01, 
        3.52810E-01, 3.63550E-01, 3.74620E-01, 3.86030E-01, 3.97790E-01, 4.09900E-01, 
        4.22380E-01, 4.35250E-01, 4.48500E-01, 4.62160E-01, 4.76240E-01, 4.90740E-01, 
        5.05690E-01, 5.21090E-01, 5.36960E-01, 5.53310E-01, 5.70160E-01, 5.87520E-01, 
        6.05420E-01, 6.23850E-01, 6.42850E-01, 6.62430E-01, 6.82600E-01, 7.03390E-01, 
        7.24810E-01, 7.46890E-01, 7.69630E-01, 7.93070E-01, 8.17230E-01, 8.42110E-01, 
        8.67760E-01, 8.94190E-01, 9.21420E-01, 9.49480E-01, 9.78400E-01, 1.00820E+00, 
        1.03890E+00, 1.07050E+00, 1.10310E+00, 1.13670E+00, 1.17140E+00, 1.20700E+00, 
        1.24380E+00, 1.28170E+00, 1.32070E+00, 1.36090E+00, 1.40240E+00, 1.44510E+00, 
        1.48910E+00, 1.53440E+00, 1.58120E+00, 1.62930E+00, 1.67890E+00, 1.73010E+00, 
        1.78280E+00, 1.83710E+00, 1.89300E+00, 1.95060E+00, 2.01010E+00, 2.07130E+00, 
        2.13430E+00, 2.19930E+00, 2.26630E+00, 2.33530E+00, 2.40650E+00, 2.47980E+00, 
        2.55530E+00, 2.63310E+00, 2.71330E+00, 2.79590E+00, 2.88110E+00, 2.96880E+00, 
        3.05920E+00, 3.15240E+00, 3.24840E+00, 3.34730E+00, 3.44930E+00, 3.55430E+00, 
        3.66250E+00, 3.77410E+00, 3.88900E+00, 4.00750E+00, 4.12950E+00, 4.25530E+00, 
        4.38490E+00, 4.51840E+00, 4.65600E+00, 4.79780E+00, 4.94390E+00, 5.09450E+00, 
        5.24960E+00, 5.40950E+00, 5.57420E+00, 5.74400E+00, 5.91890E+00, 6.09920E+00, 
        6.28490E+00, 6.47640E+00, 6.67360E+00, 6.87680E+00, 7.08630E+00, 7.30210E+00, 
        7.52440E+00, 7.75360E+00, 7.98970E+00, 8.23310E+00, 8.48380E+00, 8.74220E+00, 
        9.00840E+00, 9.28270E+00, 9.56540E+00, 9.85670E+00, 1.01570E+01, 1.04660E+01, 
        1.07850E+01, 1.11130E+01, 1.14520E+01, 1.18010E+01, 1.21600E+01, 1.25300E+01, 
        1.29120E+01, 1.33050E+01, 1.37100E+01, 1.41280E+01, 1.45580E+01, 1.50020E+01, 
        1.54580E+01, 1.59290E+01, 1.64140E+01, 1.69140E+01, 1.74290E+01, 1.79600E+01, 
        1.85070E+01, 1.90710E+01, 1.96520E+01, 2.02500E+01
    ])*1e6

    photon_energies = np.array([
        3.00000E-01, 4.00000E-01, 5.00000E-01, 6.00000E-01, 7.00000E-01, 8.00000E-01,
        9.00000E-01, 1.00000E+00, 1.10000E+00, 1.20000E+00, 1.30000E+00, 1.40000E+00,
        1.50000E+00, 1.60000E+00, 1.70000E+00, 1.80000E+00, 1.90000E+00, 2.00000E+00,
        2.10000E+00, 2.20000E+00, 2.30000E+00, 2.40000E+00, 2.50000E+00, 2.60000E+00,
        2.70000E+00, 2.80000E+00, 2.90000E+00, 3.00000E+00, 3.10000E+00, 3.20000E+00,
        3.30000E+00, 3.40000E+00, 3.50000E+00, 3.60000E+00, 3.70000E+00, 3.80000E+00,
        3.90000E+00, 4.00000E+00, 4.10000E+00, 4.20000E+00, 4.30000E+00, 4.40000E+00,
        4.50000E+00, 4.60000E+00, 4.70000E+00, 4.80000E+00, 4.90000E+00, 5.00000E+00,
        5.50000E+00, 6.00000E+00, 6.50000E+00, 7.00000E+00, 7.50000E+00, 8.00000E+00,
        8.50000E+00, 9.00000E+00, 9.50000E+00, 1.00000E+01, 1.05000E+01, 1.10000E+01
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