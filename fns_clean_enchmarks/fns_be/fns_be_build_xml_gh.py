import openmc
import numpy as np

##############################################################################
# Analysis of Beryllium Experiment

###############################################################################
# Materials

mat1 = openmc.Material(material_id=1,   name='air')
mat1.set_density('atom/b-cm',   4.921e-05)
mat1.add_nuclide('N14', 0.7886608412924202,   'ao')
mat1.add_nuclide('O16', 0.21133915870757974,  'ao')

mat2 = openmc.Material(material_id=2, name='Beryllium')
mat2.set_density('atom/b-cm',   0.122149 )
mat2.add_nuclide('Be9',  0.9927289710133025  ,    'ao')
mat2.add_nuclide('C12',  0.0006299237839521456 ,  'ao')
mat2.add_nuclide('O16',  0.0040693555162183695 ,  'ao')
mat2.add_nuclide('Al27', 0.0023701485875583395 ,  'ao')
mat2.add_nuclide('Fe54', 1.1783584234715633e-05,  'ao')
mat2.add_nuclide('Fe56', 0.00018497707234766441,  'ao')
mat2.add_nuclide('Fe57', 4.271927287144984e-06 ,  'ao')
mat2.add_nuclide('Fe58', 5.685150990914985e-07 ,  'ao')
# print(mat2.get_mass_density())

mat3 = openmc.Material(material_id=3, name='drawer-tube')
mat3.set_density('atom/b-cm',   0.0200898)
mat3.add_nuclide('Mn55', 0.00997121925432881,   'ao')
mat3.add_nuclide('Cr50', 0.008470987544935695,  'ao')
mat3.add_nuclide('Cr52', 0.1633545628084274,    'ao')
mat3.add_nuclide('Cr53', 0.018523096125301272,  'ao')
mat3.add_nuclide('Cr54', 0.004610790689015631,  'ao')
mat3.add_nuclide('Fe54', 0.041430336359409875,  'ao')
mat3.add_nuclide('Fe56', 0.6503676787547126,    'ao')
mat3.add_nuclide('Fe57', 0.015019825961606423,  'ao')
mat3.add_nuclide('Fe58', 0.001998863105791888,  'ao')
mat3.add_nuclide('Ni58', 0.05871820932193519,   'ao')
mat3.add_nuclide('Ni60', 0.022618029628936442,  'ao')
mat3.add_nuclide('Ni61', 0.0009831938364803666, 'ao')
mat3.add_nuclide('Ni62', 0.0031349384315041144, 'ao')
mat3.add_nuclide('Ni64', 0.0007982681776143339, 'ao')

# Collect materials, export to XML
materials = openmc.Materials([mat1, mat2, mat3])
materials.export_to_xml()

###############################################################################
# Geometry

# Surfaces
one = openmc.ZPlane(surface_id=1,             z0=-10.0 )
two = openmc.ZPlane(surface_id=2,             z0=20.0  )
three = openmc.ZPlane(surface_id=3,           z0=65.54 )
four = openmc.ZPlane(surface_id=4,            z0=70.00 )
five = openmc.ZCylinder(surface_id=5,         r=31.5   )
six = openmc.ZCylinder(surface_id=6,          r=35.0   )
seven = openmc.XPlane(surface_id=7,           x0=-2.46 )
eight = openmc.XPlane(surface_id=8,           x0=2.46  )
nine = openmc.XPlane(surface_id=9,            x0=-2.54 )
ten = openmc.XPlane(surface_id=10,            x0=2.54  )
eleven = openmc.YPlane(surface_id=11,         y0=-2.46 )
twelve = openmc.YPlane(surface_id=12,         y0=2.46  )
thirtn = openmc.YPlane(surface_id=13,         y0=-2.54 )
fourtn = openmc.YPlane(surface_id=14,         y0=2.54  )
fiftn = openmc.ZPlane(surface_id=15,          z0=19.9  )
sixtn = openmc.ZPlane(surface_id=16,          z0=21.0  )
sevent = openmc.ZPlane(surface_id=17,         z0=22.0  )
eightn = openmc.ZPlane(surface_id=18,         z0=23.4  )
nintn = openmc.ZPlane(surface_id=19,          z0=24.0  )
twenty = openmc.ZPlane(surface_id=20,         z0=25.0  )
twentyone = openmc.ZPlane(surface_id=21,      z0=27.0  )
twentytwo = openmc.ZPlane(surface_id=22,      z0=29.0  )
twentythree = openmc.ZPlane(surface_id=23,    z0=31.0  )
twentyfour = openmc.ZPlane(surface_id=24,     z0=33.0  )
twentyfive = openmc.ZPlane(surface_id=25,     z0=35.0  )
twentysix = openmc.ZPlane(surface_id=26,      z0=37.0  )
twentyseven = openmc.ZPlane(surface_id=27,    z0=39.0  )
twentyeight = openmc.ZPlane(surface_id=28,    z0=41.0  )
twentynine = openmc.ZPlane(surface_id=29,     z0=43.0  )
thirty = openmc.ZPlane(surface_id=30,         z0=45.0  )
thirtyone = openmc.ZPlane(surface_id=31,      z0=48.0  )
thirtytwo = openmc.ZPlane(surface_id=32,      z0=51.0  )
thirtythree = openmc.ZPlane(surface_id=33,    z0=54.0  )
thirtyfour = openmc.ZPlane(surface_id=34,     z0=57.0  )
thirtyfive = openmc.ZPlane(surface_id=35,     z0=60.0  )
thirtysix = openmc.ZPlane(surface_id=36,      z0=63.0  )
thirtyseven = openmc.ZPlane(surface_id=37,    z0=66.0  )
thirtyeight = openmc.Sphere(surface_id=38,    x0=0.0, y0=0.0, z0=0.0, r=0.5)
thirtynine = openmc.ZPlane(surface_id=39,     z0=-15.0 )
forty = openmc.ZPlane(surface_id=40,          z0=75.0  )

one.boundary_type = 'vacuum'
forty.boundary_type = 'vacuum'
six.boundary_type = 'vacuum'

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
cell11 = openmc.Cell(cell_id=11, name='cell 11')
cell12 = openmc.Cell(cell_id=12, name='cell 12')
cell13 = openmc.Cell(cell_id=13, name='cell 13')
cell14 = openmc.Cell(cell_id=14, name='cell 14')
cell15 = openmc.Cell(cell_id=15, name='cell 15')
cell16 = openmc.Cell(cell_id=16, name='cell 16')
cell17 = openmc.Cell(cell_id=17, name='cell 17')
cell18 = openmc.Cell(cell_id=18, name='cell 18')
cell19 = openmc.Cell(cell_id=19, name='cell 19')
cell20 = openmc.Cell(cell_id=20, name='cell 20')
cell21 = openmc.Cell(cell_id=21, name='cell 21')
cell22 = openmc.Cell(cell_id=22, name='cell 22')
cell23 = openmc.Cell(cell_id=23, name='cell 23') 
cell24 = openmc.Cell(cell_id=24, name='cell 24')
cell25 = openmc.Cell(cell_id=25, name='cell 25')
cell26 = openmc.Cell(cell_id=26, name='cell 26')
cell27 = openmc.Cell(cell_id=27, name='cell 27')
cell28 = openmc.Cell(cell_id=28, name='cell 28')
cell29 = openmc.Cell(cell_id=29, name='cell 29')
cell30 = openmc.Cell(cell_id=30, name='cell 30')

# Regions
cell1.region = +one & -two & -five & +thirtyeight & ~ (+fiftn & -two & +seven & -eight & +eleven & -twelve)
cell2.region = +three & -four & -five & ~ (+three & -thirtyseven & +seven & -eight & +eleven & -twelve)
cell3.region = (+two & -three & -five & -nine) | (+two & -three & -five & +ten) | (+two & -three & -five & -thirtn) | (+two & -three & -five & +fourtn)
cell4.region = (+two & -three & +nine & -seven & +thirtn & -fourtn) | (+two & -three & +eight & -ten & +thirtn & -fourtn) | (+two & -three & +seven & -eight & +thirtn & -eleven) | (+two & -three & +seven & -eight & +twelve & -fourtn) 
cell5.region =   +fiftn       & -two            & +seven & -eight & +eleven & -twelve
cell6.region =   +two         & -sixtn          & +seven & -eight & +eleven & -twelve
cell7.region =   +sixtn       & -sevent         & +seven & -eight & +eleven & -twelve
cell8.region =   +sevent      & -eightn         & +seven & -eight & +eleven & -twelve
cell9.region =   +eightn      & -nintn          & +seven & -eight & +eleven & -twelve
cell10.region =  +nintn       & -twenty         & +seven & -eight & +eleven & -twelve
cell11.region =  +twenty      & -twentyone      & +seven & -eight & +eleven & -twelve
cell12.region =  +twentyone   & -twentytwo      & +seven & -eight & +eleven & -twelve
cell13.region =  +twentytwo   & -twentythree    & +seven & -eight & +eleven & -twelve
cell14.region =  +twentythree & -twentyfour     & +seven & -eight & +eleven & -twelve
cell15.region =  +twentyfour  & -twentyfive     & +seven & -eight & +eleven & -twelve
cell16.region =  +twentyfive  & -twentysix      & +seven & -eight & +eleven & -twelve
cell17.region =  +twentysix   & -twentyseven    & +seven & -eight & +eleven & -twelve
cell18.region =  +twentyseven & -twentyeight    & +seven & -eight & +eleven & -twelve
cell19.region =  +twentyeight & -twentynine     & +seven & -eight & +eleven & -twelve
cell20.region =  +twentynine  & -thirty         & +seven & -eight & +eleven & -twelve
cell21.region =  +thirty      & -thirtyone      & +seven & -eight & +eleven & -twelve
cell22.region =  +thirtyone   & -thirtytwo      & +seven & -eight & +eleven & -twelve
cell23.region =  +thirtytwo   & -thirtythree    & +seven & -eight & +eleven & -twelve
cell24.region =  +thirtythree & -thirtyfour     & +seven & -eight & +eleven & -twelve
cell25.region =  +thirtyfour  & -thirtyfive     & +seven & -eight & +eleven & -twelve
cell26.region =  +thirtyfive  & -thirtysix      & +seven & -eight & +eleven & -twelve
cell27.region =  +thirtysix   & -three          & +seven & -eight & +eleven & -twelve
cell28.region =  +three       & -thirtyseven    & +seven & -eight & +eleven & -twelve 
cell29.region =  -thirtyeight
cell30.region =  (+thirtynine & -one & -six) | (+one & -four & -six & +five) | (+four & -forty & -six)

# cell materials
cell1.fill = mat1
cell2.fill = mat1
cell3.fill = mat2
cell4.fill = mat3
cell5.fill = mat1
cell6.fill = mat2
cell7.fill = mat2
cell8.fill = mat2
cell9.fill = mat2
cell10.fill = mat2
cell11.fill = mat2
cell12.fill = mat2
cell13.fill = mat2
cell14.fill = mat2
cell15.fill = mat2
cell16.fill = mat2
cell17.fill = mat2
cell18.fill = mat2
cell19.fill = mat2
cell20.fill = mat2
cell21.fill = mat2
cell22.fill = mat2
cell23.fill = mat2
cell24.fill = mat2
cell25.fill = mat2
cell26.fill = mat2
cell27.fill = mat2
cell28.fill = mat1
cell29.fill = mat1

# Create a geometry and export to XML
geometry = openmc.Geometry([cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9,cell10,cell11,cell12,cell13,cell14,cell15,cell16,cell17,cell18,cell19,cell20,cell21,cell22,cell23,cell24,cell25,cell26,cell27,cell28,cell29,cell30])
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Indicate how many particles to run
settings = openmc.Settings()
settings.run_mode = 'fixed source'
settings.batches = 10000000
settings.particles = 200
settings.survival_biasing = True 
settings.photon_transport = True

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
source.space = openmc.stats.Point(xyz=(0.0,0.0,0.0))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Tabular(x,p,interpolation='histogram')
settings.source = openmc.source.Source()
settings.export_to_xml()

###############################################################################
# Tallies gamma heating   
cell_filter = openmc.CellFilter([cell5,cell6,cell7,cell8,cell9,cell10,cell11,cell12,cell13,cell14,cell15,cell16,cell17,cell18,cell19,cell20,cell21,cell22,cell23,cell24,cell25,cell26])
particle_p = openmc.ParticleFilter(['neutron','photon','electron','positron'])
tally_p = openmc.Tally(name="Gamma heating")  
tally_p.filters = [cell_filter,particle_p]
tally_p.scores = ['heating']

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([tally_p])
tallies.export_to_xml()

openmc.run()