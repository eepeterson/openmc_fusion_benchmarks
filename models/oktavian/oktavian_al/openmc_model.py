#%%
import openmc
import numpy as np

#%%

# MATERIALS

# 
mat_1 = openmc.Material(material_id=1)
mat_1.set_density('g/cm3', 1.223)
mat_1.add_nuclide('Al27', 0.9975488, 'ao')
mat_1.add_nuclide('Si28', 0.001329808, 'ao')
mat_1.add_nuclide('Si29', 6.752131e-05, 'ao')
mat_1.add_nuclide('Si30', 4.450956e-05, 'ao')
mat_1.add_nuclide('Fe54', 5.651123e-05, 'ao')
mat_1.add_nuclide('Fe56', 0.0008871055, 'ao')
mat_1.add_nuclide('Fe57', 2.048713e-05, 'ao')
mat_1.add_nuclide('Fe58', 2.726461e-06, 'ao')
mat_1.add_nuclide('Cu63', 2.938581e-05, 'ao')
mat_1.add_nuclide('Cu65', 1.309765e-05, 'ao')
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

materials.export_to_xml()

#%%

# GEOMETRY

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
geometry = openmc.Geometry(universe)
geometry.export_to_xml()

#%%

# SETTINGS


# 
# Instantiate a Settings object
settings = openmc.Settings()
settings.batches = 100
settings.particles = 1000000000
settings.run_mode = 'fixed source'
settings.photon_transport = True  # This line is required to switch on photons tracking, other wise the photons created by the neutrons are not tracked

# Create a point source
energies = np.array([
    0.14000, 0.15079, 0.16665, 0.18418, 0.20355, 0.22496, 0.24862, 0.27476, 
    0.30366, 0.33560, 0.37089, 0.40990, 0.45301, 0.50065, 0.55331, 0.61150, 
    0.67581, 0.74689, 0.82544, 0.91225, 1.00820, 1.11420, 1.23140, 1.36090, 
    1.50400, 1.66220, 1.83700, 2.03020, 2.24380, 2.47970, 2.74050, 3.02880, 
    3.34730, 3.69930, 4.08840, 4.51840, 4.99360, 5.51880, 6.09920, 6.74060, 
    7.44960, 8.23300, 9.09890, 10.05600, 11.11300, 12.28200, 13.57400, 15.00200, 
    16.57900, 18.32300, 20.25000
])*1e6

weights = np.array([
    0.000193606, 0.000142111, 0.000212238, 0.00054769, 0.00027653, 
    0.000374911, 0.000670615, 0.000617066, 0.000660875, 0.000784287, 
    0.000956674, 0.000869708, 0.0013789, 0.00146416, 0.00170894, 0.00181295, 
    0.00192411, 0.00214897, 0.00228494, 0.00241619, 0.00272518, 0.00290644, 
    0.00279685, 0.00301846, 0.00304442, 0.00297376, 0.00295501, 0.00313515, 
    0.00288504, 0.0034731, 0.00324985, 0.00269866, 0.00254773, 0.00243325, 
    0.00221309, 0.00204398, 0.00183606, 0.00179295, 0.00173385, 0.00157818, 
    0.00156385, 0.00179906, 0.00235438, 0.00308434, 0.00508898, 0.0135129, 
    0.652778, 0.248384, 0.000320765, 8.02381e-05, 0.0
])

source = openmc.IndependentSource()
source.space = openmc.stats.Point((0, 0, 0))
source.angle = openmc.stats.Isotropic()
source.particle = 'neutron' 
source.energy = openmc.stats.Tabular(energies, weights, interpolation='histogram')
settings.source = source

settings.export_to_xml()
#%%

# TALLIES

# 
# setup the filters for the cell tally
vessel_surface_filter = openmc.SurfaceFilter(surf_6) # outer radius of outer steel vessel
neutron_particle_filter = openmc.ParticleFilter(['neutron'])
photon_particle_filter = openmc.ParticleFilter(['photon'])

neutron_energies = np.array([
    0.09712, 0.10109, 0.10521, 0.10950,
    0.11397, 0.11862, 0.12347, 0.12850, 0.13375, 0.13921,
    0.14489, 0.15080, 0.15696, 0.16336, 0.17003, 0.17697,
    0.18419, 0.19171, 0.19953, 0.20767, 0.21615, 0.22497,
    0.23415, 0.24371, 0.25365, 0.26400, 0.27478, 0.28599,
    0.29766, 0.30981, 0.32245, 0.33561, 0.34931, 0.36357,
    0.37840, 0.39385, 0.40992, 0.42665, 0.44406, 0.46218,
    0.48105, 0.50068, 0.52111, 0.54238, 0.56451, 0.58755,
    0.61153, 0.63648, 0.66246, 0.68950, 0.71763, 0.74692,
    0.77740, 0.80913, 0.84215, 0.87652, 0.91229, 0.94952,
    0.98827, 1.02860, 1.07060, 1.11430, 1.15980, 1.20710,
    1.25630, 1.30760, 1.36100, 1.41650, 1.47430, 1.53450,
    1.59710, 1.66230, 1.73010, 1.80080, 1.87420, 1.95070,
    2.03030, 2.11320, 2.19940, 2.28920, 2.38260, 2.47990,
    2.58110, 2.68640, 2.79600, 2.91010, 3.02890, 3.15250,
    3.28120, 3.41510, 3.55450, 3.69950, 3.85050, 4.00760,
    4.17120, 4.34140, 4.51860, 4.70300, 4.89490, 5.09470,
    5.30260, 5.51900, 5.74430, 5.97870, 6.22270, 6.47660,
    6.74100, 7.01610, 7.30240, 7.60040, 7.91060, 8.23340,
    8.56940, 8.91920, 9.28320, 9.66200, 10.05600, 10.46700,
    10.89400, 11.33900, 11.80100, 12.28300, 12.78400, 13.30600,
    13.84900, 14.41400, 15.00200, 15.61500, 16.25200, 16.91500,
    17.60500, 18.32400, 19.07200, 19.85000, 20.66000
])*1e6

photon_energies = np.array([
    0.50000, 0.60000, 0.70000, 0.80000,
    0.90000, 1.00000, 1.10000, 1.20000, 1.30000, 1.40000,
    1.50000, 1.60000, 1.70000, 1.80000, 1.90000, 2.00000,
    2.10000, 2.20000, 2.30000, 2.40000, 2.50000, 2.60000,
    2.70000, 2.80000, 2.90000, 3.00000, 3.10000, 3.20000,
    3.30000, 3.40000, 3.50000, 3.60000, 3.70000, 3.80000,
    3.90000, 4.00000, 4.10000, 4.20000, 4.30000, 4.40000,
    4.50000, 4.60000, 4.70000, 4.80000, 4.90000, 5.00000,
    5.50000, 6.00000, 6.50000, 7.00000, 7.50000, 8.00000,
    8.50000, 9.00000, 9.50000, 10.00000, 10.50000,
    11.00000
])*1e6

energy_neutron_filter = openmc.EnergyFilter(neutron_energies)
energy_photon_filter = openmc.EnergyFilter(photon_energies)

#creates an empty tally object
tallies = openmc.Tallies()

# create the tally
neutron_surface_spectra_tally = openmc.Tally(name='vessel_neutron_spectra_tally')
neutron_surface_spectra_tally.scores = ['current']
neutron_surface_spectra_tally.filters = [vessel_surface_filter, neutron_particle_filter, energy_neutron_filter]
tallies.append(neutron_surface_spectra_tally)

photon_surface_spectra_tally = openmc.Tally(name='vessel_photon_spectra_tally')
photon_surface_spectra_tally.scores = ['current']
photon_surface_spectra_tally.filters = [vessel_surface_filter, photon_particle_filter, energy_photon_filter]
tallies.append(photon_surface_spectra_tally)

tallies.export_to_xml()