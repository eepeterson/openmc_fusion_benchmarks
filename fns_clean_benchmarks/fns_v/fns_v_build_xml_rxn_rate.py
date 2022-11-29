import openmc
import numpy as np

##############################################################################
# Analysis of Vanadium Experiment

###############################################################################
# Materials

mat1 = openmc.Material(material_id=1, name='Vanadium')
mat1.set_density('g/cc', 6.033       )
mat1.add_element('V', 0.99777,   'wo')
mat1.add_element('Al', 0.073e-2, 'wo')
mat1.add_element('Si', 0.108e-2, 'wo')
mat1.add_element('Fe', 0.042e-2, 'wo')
mat1.add_element('Nb', 1e-8,     'wo')
mat1.add_element('In', 1e-8,     'wo')
mat1.add_element('Au', 1e-8,     'wo')
# print(mat1.get_nuclide_atom_densities())

mat2 = openmc.Material(material_id=2, name='Graphite')
mat2.set_density('g/cc', 1.625     )
mat2.add_nuclide('C12', 1.0,   'ao')

mat3 = openmc.Material(material_id=3, name='air')
mat3.set_density('atom/b-cm',   4.921e-05)
mat3.add_nuclide('N14',  0.231, 'wo')
mat3.add_nuclide('O16',  0.756, 'wo')
mat3.add_nuclide('Ar40', 0.013, 'wo')
mat3.add_element('Al', 1e-8,    'wo')
mat3.add_element('Nb', 1e-8,    'wo')
mat3.add_element('In', 1e-8,    'wo')
mat3.add_element('Au', 1e-8,    'wo')
# print(mat3.get_nuclide_atom_densities())

# Collect the materials together and export to XML
materials = openmc.Materials([mat1, mat2, mat3])
materials.export_to_xml()

###############################################################################
# Geometry

# Surfaces
one = openmc.ZPlane(surface_id=1,            z0=-21.0  )
two = openmc.ZPlane(surface_id=2,            z0=-2.0   )
three = openmc.ZPlane(surface_id=3,          z0=-0.10  )
four = openmc.ZPlane(surface_id=4,           z0=0.00   )
five = openmc.ZPlane(surface_id=5,           z0=7.37   )
six = openmc.ZPlane(surface_id=6,            z0=7.87   )
seven = openmc.ZPlane(surface_id=7,          z0=10.00  )
eight = openmc.ZPlane(surface_id=8,          z0=17.53  )
nine = openmc.ZPlane(surface_id=9,           z0=18.03  )
ten = openmc.ZPlane(surface_id=10,           z0=25.40  )
eleve = openmc.ZPlane(surface_id=11,         z0=30.48  )
twentyone = openmc.XPlane(surface_id=21,     x0=-12.70 )
twentytwo = openmc.XPlane(surface_id=22,     x0=12.70  )
twentythree = openmc.YPlane(surface_id=23,   y0=-12.70 )
twentyfour = openmc.YPlane(surface_id=24,    y0=12.70  )
twentyfive = openmc.XPlane(surface_id=25,    x0=-17.78 )
twentysix = openmc.XPlane(surface_id=26,     x0=17.78  )
twentyseven = openmc.YPlane(surface_id=27,   y0=-17.78 )
twentyeight = openmc.YPlane(surface_id=28,   y0=17.78  )
twentynine = openmc.ZCylinder(surface_id=29, r=1.00    )

one.boundary_type = 'vacuum'
eleve.boundary_type = 'vacuum'
twentyfive.boundary_type = 'vacuum'
twentysix.boundary_type = 'vacuum'
twentyseven.boundary_type = 'vacuum' 
twentyeight.boundary_type = 'vacuum' 

# Cells
cell1 = openmc.Cell(cell_id=1, name='Cell 1')  
cell2 = openmc.Cell(cell_id=2, name='cell 2') 
cell3 = openmc.Cell(cell_id=3, name='cell 3')
cell4 = openmc.Cell(cell_id=4, name='cell 4')
cell5 = openmc.Cell(cell_id=5, name='cell 5')
cell6 = openmc.Cell(cell_id=6, name='cell 6')
cell7 = openmc.Cell(cell_id=7, name='cell 7')
cell8 = openmc.Cell(cell_id=8, name='cell 8')
cell9 = openmc.Cell(cell_id=9, name='cell 9')
cell10 = openmc.Cell(cell_id=10, name='cell 10')

# cell regions
cell1.region = +twentyfive & -twentysix & +twentyseven & -twentyeight & +one & -four & ~ (+two & -four & -twentynine) 
cell2.region = +two & -three & -twentynine
cell3.region = +three & -four & -twentynine
cell4.region = +five & -six & -twentynine
cell5.region = +eight & -nine & -twentynine
cell6.region = +twentyone & -twentytwo  & +twentythree & -twentyfour  & +four  & -seven & ~ (+five & -six &  -twentynine)
cell7.region = +twentyone & -twentytwo  & +twentythree & -twentyfour  & +seven & -ten & ~ (+eight & -nine & -twentynine)
cell8.region = +twentyfive & -twentysix & +twentyseven & -twentyeight & +four  & -seven & ~ (+twentyone & -twentytwo & +twentythree & -twentyfour & +four & -seven)
cell9.region = +twentyfive & -twentysix & +twentyseven & -twentyeight & +seven & -eleve & ~ (+twentyone & -twentytwo & +twentythree & -twentyfour & +seven & -ten)
cell10.region =  -twentyfive | +twentysix | -twentyseven | +twentyeight | -one | +eleve 

# cell materials
cell1.fill = mat3
cell2.fill = mat3
cell3.fill = mat3
cell4.fill = mat1
cell5.fill = mat1
cell6.fill = mat1
cell7.fill = mat1
cell8.fill = mat2
cell9.fill = mat2

# export geometry to XML
geometry = openmc.Geometry([cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9,cell10])
geometry.export_to_xml()

###############################################################################
# Settings

# Number of particles to run
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.batches = 500000
settings.particles = 200
settings.survival_biasing = True
# settings.photon_transport = True 

# Create source distribution 
# upper neutron energy bins
x = 1e6*np.array([1.0010e-11,  3.2241e-07,
    5.3156e-07,  8.7640e-07,  1.4449e-06,  2.3823e-06,  3.9278e-06,
    6.4758e-06,  1.0677e-05,  1.7603e-05,  2.9023e-05,  4.7850e-05,
    7.8891e-05,  1.3007e-04,  2.1445e-04,  3.5357e-04,  5.8293e-04,
    9.6110e-04,  1.2341e-03,  1.5846e-03,  2.0346e-03,  2.6125e-03,
    3.3546e-03,  4.3073e-03,  5.5307e-03,  7.1016e-03,  9.1186e-03,
    1.1709e-02,  1.5034e-02,  1.9304e-02,  2.1874e-02,  2.4787e-02,
    2.8087e-02,  3.1827e-02,  3.6065e-02,  4.0867e-02,  4.6308e-02,
    5.2474e-02,  5.9461e-02,  6.7378e-02,  7.6349e-02,  8.6515e-02,
    9.8035e-02,  1.1109e-01,  1.2588e-01,  1.4264e-01,  1.6163e-01,
    1.8315e-01,  2.0754e-01,  2.3517e-01,  2.6649e-01,  3.0197e-01,
    3.4217e-01,  3.8774e-01,  4.3936e-01,  4.9786e-01,  5.6415e-01,
    6.3927e-01,  7.2438e-01,  8.2084e-01,  9.3013e-01,  1.0540e+00,
    1.1943e+00,  1.3533e+00,  1.5335e+00,  1.7377e+00,  1.8498e+00,
    1.9691e+00,  2.0961e+00,  2.2313e+00,  2.3752e+00,  2.5284e+00,
    2.6914e+00,  2.8650e+00,  3.0498e+00,  3.2465e+00,  3.4559e+00,
    3.6787e+00,  3.9160e+00,  4.1686e+00,  4.4374e+00,  4.7236e+00,
    5.0282e+00,  5.3525e+00,  5.6978e+00,  6.0652e+00,  6.4564e+00,
    6.8728e+00,  7.3161e+00,  7.7879e+00,  8.2902e+00,  8.8249e+00,
    9.3940e+00,  9.9999e+00,  1.0157e+01,  1.0317e+01,  1.0480e+01,
    1.0645e+01,  1.0812e+01,  1.0983e+01,  1.1156e+01,  1.1331e+01,
    1.1510e+01,  1.1691e+01,  1.1875e+01,  1.2062e+01,  1.2252e+01,
    1.2445e+01,  1.2641e+01,  1.2840e+01,  1.3042e+01,  1.3248e+01,
    1.3456e+01,  1.3668e+01,  1.3883e+01,  1.4102e+01,  1.4324e+01,
    1.4550e+01,  1.4779e+01,  1.5012e+01,  1.5248e+01,  1.5488e+01])

# bin probabilities       
p = np.array([0.0, 1.5142e-07,
    2.2732e-09,  4.2225e-09,  7.4848e-09,  1.4264e-08,  8.3975e-08,
    1.8398e-07,  2.2450e-07,  1.3922e-07,  1.6817e-07,  2.9754e-07,
    3.8068e-06,  3.0541e-06,  2.2612e-06,  6.9372e-06,  7.2049e-06,
    8.7622e-06,  7.8013e-06,  1.4320e-05,  1.1820e-05,  1.6544e-05,
    1.4791e-05,  1.7624e-05,  2.8404e-05,  2.4899e-05,  3.7633e-05,
    4.4237e-05,  4.6320e-05,  6.1572e-05,  3.7185e-05,  5.3362e-05,
    4.8831e-05,  5.0292e-05,  5.7202e-05,  6.9230e-05,  8.0602e-05,
    8.3190e-05,  9.7450e-05,  1.0531e-04,  1.2632e-04,  1.4874e-04,
    1.7906e-04,  3.7225e-04,  4.9933e-04,  5.3824e-04,  6.0762e-04,
    7.0593e-04,  8.0965e-04,  9.5392e-04,  1.0785e-03,  1.2232e-03,
    1.3867e-03,  1.5803e-03,  1.6473e-03,  1.8238e-03,  2.0605e-03,
    2.2042e-03,  2.3040e-03,  2.5211e-03,  2.5709e-03,  2.5872e-03,
    2.5765e-03,  2.7699e-03,  2.8528e-03,  2.5945e-03,  1.3898e-03,
    1.4298e-03,  1.3270e-03,  1.3489e-03,  1.3820e-03,  1.4312e-03,
    1.3760e-03,  1.4329e-03,  1.4558e-03,  1.3518e-03,  1.4053e-03,
    1.2861e-03,  1.2741e-03,  1.1711e-03,  1.1937e-03,  1.0563e-03,
    1.0018e-03,  8.8451e-04,  7.9827e-04,  7.9293e-04,  7.5872e-04,
    6.9228e-04,  6.2956e-04,  5.1710e-04,  5.0750e-04,  5.1007e-04,
    4.1280e-04,  3.5649e-04,  9.0768e-05,  8.2287e-05,  9.2862e-05,
    9.1407e-05,  9.3708e-05,  7.9567e-05,  8.8737e-05,  8.7841e-05,
    1.1227e-04,  1.6798e-04,  1.5985e-04,  1.6563e-04,  2.1025e-04,
    4.1363e-04,  7.4899e-04,  7.8183e-04,  5.1771e-04,  4.5938e-04,
    4.6458e-04,  9.1020e-04,  2.6083e-03,  9.5007e-04,  5.1474e-03,
    3.0897e-02,  2.3565e-01,  4.0901e-01,  2.2296e-01,  1.4419e-01])

source = openmc.Source
source.particle = 'neutron'
source.space = openmc.stats.Point(xyz=(0.0,0.0,-20.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Tabular(x,p,interpolation='histogram')
settings.source = openmc.source.Source()
settings.export_to_xml()

###############################################################################
# Tallies

# neutron energy spectrum tally
cell_filter = openmc.CellFilter([cell3, cell4, cell5])
particle_n = openmc.ParticleFilter(['neutron'])
tally_1 = openmc.Tally(name="Al-27 (n,alpha)")
tally_1.nuclides = ['Al27']
tally_1.filters = [cell_filter,particle_n]
tally_1.scores = ['(n,a)']

tally_2 = openmc.Tally(name="Nb-93 (n,2n)")
tally_2.nuclides = ['Nb93']
tally_2.filters = [cell_filter,particle_n]
tally_2.scores = ['(n,2n)']

tally_3 = openmc.Tally(name="In115 (n,n')")
tally_3.nuclides = ['In115']
tally_3.filters = [cell_filter,particle_n]
tally_3.scores = ['(n,nc)']

tally_4 = openmc.Tally(name="Au-197 (n,gamma)")
tally_4.nuclides = ['Au197']
tally_4.filters = [cell_filter,particle_n]
tally_4.scores = ['(n,gamma)']

# Collect Tallies, export to XML
tallies = openmc.Tallies([tally_1,tally_2,tally_3,tally_4])
tallies.export_to_xml()

openmc.run()