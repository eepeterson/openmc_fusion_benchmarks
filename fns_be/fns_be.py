#!/usr/bin/env python3
import numpy as np
import argparse

import openmc


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--batches', type=int, default=100)
    parser.add_argument('-p', '--particles', type=int, default=10000)
    parser.add_argument('-s', '--threads', type=int)
    parser.add_argument('-c', '--cwd', type=str)
    group = parser.add_argument_group('tallies')
    group.add_argument('-n', '--neutron_spectrum', action='store_true',
                       default=False)
    group.add_argument('-g', '--gamma_heating', action='store_true',
                       default=False)
    group.add_argument('-r', '--reaction_rates', action='store_true',
                       default=False)

    args = parser.parse_args()

    # If neither gamma heating or rxn rate tallies are desired, default to only
    # the neutron spectrum tallies
    if not (args.gamma_heating or args.reaction_rates):
        args.neutron_spectrum = True

    return args


def main():
    """Analysis of Beryllium Experiment"""

    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    # Build materials
    air = openmc.Material(name='air')
    air.add_nuclide('N14', 0.7886608412924202, 'ao')
    air.add_nuclide('O16', 0.21133915870757974, 'ao')

    beryllium = openmc.Material(name='Beryllium')
    beryllium.add_nuclide('Be9', 0.9927289710133025, 'ao')
    beryllium.add_nuclide('C12', 0.0006299237839521456, 'ao')
    beryllium.add_nuclide('O16', 0.0040693555162183695, 'ao')
    beryllium.add_nuclide('Al27', 0.0023701485875583395, 'ao')
    beryllium.add_nuclide('Fe54', 1.1783584234715633e-05, 'ao')
    beryllium.add_nuclide('Fe56', 0.00018497707234766441, 'ao')
    beryllium.add_nuclide('Fe57', 4.271927287144984e-06, 'ao')
    beryllium.add_nuclide('Fe58', 5.685150990914985e-07, 'ao')

    tube = openmc.Material(name='drawer-tube')
    tube.add_nuclide('Mn55', 0.00997121925432881, 'ao')
    tube.add_nuclide('Cr50', 0.008470987544935695, 'ao')
    tube.add_nuclide('Cr52', 0.1633545628084274, 'ao')
    tube.add_nuclide('Cr53', 0.018523096125301272, 'ao')
    tube.add_nuclide('Cr54', 0.004610790689015631, 'ao')
    tube.add_nuclide('Fe54', 0.041430336359409875, 'ao')
    tube.add_nuclide('Fe56', 0.6503676787547126, 'ao')
    tube.add_nuclide('Fe57', 0.015019825961606423, 'ao')
    tube.add_nuclide('Fe58', 0.001998863105791888, 'ao')
    tube.add_nuclide('Ni58', 0.05871820932193519, 'ao')
    tube.add_nuclide('Ni60', 0.022618029628936442, 'ao')
    tube.add_nuclide('Ni61', 0.0009831938364803666, 'ao')
    tube.add_nuclide('Ni62', 0.0031349384315041144, 'ao')
    tube.add_nuclide('Ni64', 0.0007982681776143339, 'ao')

    if args.reaction_rates:
        # Add reaction rate elements to materials in small quantities
        air.add_element('Al', 1e-8, 'ao')
        air.add_element('Nb', 1e-8, 'ao')
        air.add_element('In', 1e-8, 'ao')
        air.add_element('Au', 1e-8, 'ao')
        air.add_element('Fe', 1e-8, 'ao')
        air.add_element('Zr', 1e-8, 'ao')
        air.add_element('Ni', 1e-8, 'ao')
        air.add_element('Ti', 1e-8, 'ao')
        air.add_nuclide('U235', 1e-8, 'ao')
        air.add_nuclide('Li6', 1e-8, 'ao')
        air.add_nuclide('Be9', 1e-8, 'ao')
        beryllium.add_element('Nb', 1e-8, 'ao')
        beryllium.add_element('In', 1e-8, 'ao')
        beryllium.add_element('Au', 1e-8, 'ao')
        beryllium.add_element('Zr', 1e-8, 'ao')
        beryllium.add_element('Ni', 1e-8, 'ao')
        beryllium.add_element('Ti', 1e-8, 'ao')
        beryllium.add_nuclide('U235', 1e-8, 'ao')
        beryllium.add_nuclide('Li6', 1e-8, 'ao')
        tube.add_element('Al', 1e-8, 'ao')
        tube.add_element('Nb', 1e-8, 'ao')
        tube.add_element('In', 1e-8, 'ao')
        tube.add_element('Au', 1e-8, 'ao')
        tube.add_element('Zr', 1e-8, 'ao')
        tube.add_element('Ti', 1e-8, 'ao')
        tube.add_nuclide('U235', 1e-8, 'ao')
        tube.add_nuclide('Li6', 1e-8, 'ao')

    air.set_density('atom/b-cm', 4.921e-05)
    beryllium.set_density('atom/b-cm', 0.122149)
    tube.set_density('atom/b-cm', 0.0200898)

    ############################################################################
    # Build Geometry

    # Surfaces
    s1 = openmc.ZPlane(-10.0, boundary_type='vacuum')
    s2 = openmc.ZPlane(20.0)
    s3 = openmc.ZPlane(65.54)
    s4 = openmc.ZPlane(70.00)
    s5 = openmc.ZCylinder(r=31.5)
    s6 = openmc.ZCylinder(r=35.0, boundary_type='vacuum')
    s7 = openmc.XPlane(-2.46)
    s8 = openmc.XPlane(2.46)
    s9 = openmc.XPlane(-2.54)
    s10 = openmc.XPlane(2.54)
    s11 = openmc.YPlane(-2.46)
    s12 = openmc.YPlane(2.46)
    s13 = openmc.YPlane(-2.54)
    s14 = openmc.YPlane(2.54)
    s15 = openmc.ZPlane(19.9)
    s16 = openmc.ZPlane(21.0)
    s17 = openmc.ZPlane(22.0)
    s18 = openmc.ZPlane(23.4)
    s19 = openmc.ZPlane(24.0)
    s20 = openmc.ZPlane(25.0)
    s21 = openmc.ZPlane(27.0)
    s22 = openmc.ZPlane(29.0)
    s23 = openmc.ZPlane(31.0)
    s24 = openmc.ZPlane(33.0)
    s25 = openmc.ZPlane(35.0)
    s26 = openmc.ZPlane(37.0)
    s27 = openmc.ZPlane(39.0)
    s28 = openmc.ZPlane(41.0)
    s29 = openmc.ZPlane(43.0)
    s30 = openmc.ZPlane(45.0)
    s31 = openmc.ZPlane(48.0)
    s32 = openmc.ZPlane(51.0)
    s33 = openmc.ZPlane(54.0)
    s34 = openmc.ZPlane(57.0)
    s35 = openmc.ZPlane(60.0)
    s36 = openmc.ZPlane(63.0)
    s37 = openmc.ZPlane(66.0)
    s38 = openmc.Sphere(r=0.5)
    s39 = openmc.ZPlane(-15.0)
    s40 = openmc.ZPlane(75.0, boundary_type='vacuum')

    # Helper Regions
    xybox = +s7 & -s8 & +s11 & -s12
    zcyl = +s2 & -s3 & -s5
    region1 = +s1 & -s2 & -s5 & +s38 & ~(+s15 & -s2 & xybox)
    region2 = +s3 & -s4 & -s5 & ~(+s3 & -s37 & xybox)
    region3 = (zcyl & -s9) | (zcyl & +s10) | (zcyl & -s13) | (zcyl & +s14)
    region4 = (+s2 & -s3 & +s9 & -s7 & +s13 & -s14)
    region4 |= (+s2 & -s3 & +s8 & -s10 & +s13 & -s14)
    region4 |= (+s2 & -s3 & +s7 & -s8 & +s13 & -s11)
    region4 |= (+s2 & -s3 & +s7 & -s8 & +s12 & -s14)
    region30 =  (+s39 & -s1 & -s6) | (+s1 & -s4 & -s6 & +s5) | (+s4 & -s40 & -s6)

    # Cells
    cell1 = openmc.Cell(region=region1, fill=air)
    cell2 = openmc.Cell(region=region2, fill=air)
    cell3 = openmc.Cell(region=region3, fill=beryllium)
    cell4 = openmc.Cell(region=region4, fill=tube)
    cell5 = openmc.Cell(region=+s15 & -s2 & xybox, fill=air)
    cell6 = openmc.Cell(region=+s2 & -s16 & xybox, fill=beryllium)
    cell7 = openmc.Cell(region=+s16 & -s17 & xybox, fill=beryllium)
    cell8 = openmc.Cell(region=+s17 & -s18 & xybox, fill=beryllium)
    cell9 = openmc.Cell(region=+s18 & -s19 & xybox, fill=beryllium)
    cell10 = openmc.Cell(region=+s19 & -s20 & xybox, fill=beryllium)
    cell11 = openmc.Cell(region=+s20 & -s21 & xybox, fill=beryllium)
    cell12 = openmc.Cell(region=+s21 & -s22 & xybox, fill=beryllium)
    cell13 = openmc.Cell(region=+s22 & -s23 & xybox, fill=beryllium)
    cell14 = openmc.Cell(region=+s23 & -s24 & xybox, fill=beryllium)
    cell15 = openmc.Cell(region=+s24 & -s25 & xybox, fill=beryllium)
    cell16 = openmc.Cell(region=+s25 & -s26 & xybox, fill=beryllium)
    cell17 = openmc.Cell(region=+s26 & -s27 & xybox, fill=beryllium)
    cell18 = openmc.Cell(region=+s27 & -s28 & xybox, fill=beryllium)
    cell19 = openmc.Cell(region=+s28 & -s29 & xybox, fill=beryllium)
    cell20 = openmc.Cell(region=+s29 & -s30 & xybox, fill=beryllium)
    cell21 = openmc.Cell(region=+s30 & -s31 & xybox, fill=beryllium)
    cell22 = openmc.Cell(region=+s31 & -s32 & xybox, fill=beryllium)
    cell23 = openmc.Cell(region=+s32 & -s33 & xybox, fill=beryllium)
    cell24 = openmc.Cell(region=+s33 & -s34 & xybox, fill=beryllium)
    cell25 = openmc.Cell(region=+s34 & -s35 & xybox, fill=beryllium)
    cell26 = openmc.Cell(region=+s35 & -s36 & xybox, fill=beryllium)
    cell27 = openmc.Cell(region=+s36 & -s3 & xybox, fill=beryllium)
    cell28 = openmc.Cell(region=+s3 & -s37 & xybox, fill=air)
    cell29 = openmc.Cell(region=-s38, fill=air)
    cell30 = openmc.Cell(region=region30)

    cells = [cell1, cell2, cell3, cell4, cell5, cell6, cell7, cell8, cell9, cell10,
             cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18,
             cell19, cell20, cell21, cell22, cell23, cell24, cell25, cell26,
             cell27, cell28, cell29, cell30
             ]

    for c in cells:
        c.name = f'cell {c.id}'

    # Add Geometry to Model
    model.geometry = openmc.Geometry(cells)

    ###############################################################################
    # Define problem settings

    # Indicate how many particles to run
    settings = openmc.Settings(run_mode='fixed source')
    settings.batches = args.batches
    settings.particles = args.particles

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

    source = openmc.Source()
    source.particle = 'neutron'
    source.space = openmc.stats.Point()
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Tabular(x, p, interpolation='histogram')
    settings.source = source
    settings.output = {'tallies': False}

    ###############################################################################
    # Build tallies based on the desired run
    model.tallies = openmc.Tallies()
    cell_filter = openmc.CellFilter(cells[4:26])
    usrcwd = args.cwd
    cwd = ''

    if args.neutron_spectrum:
        # neutron energy spectrum tally   
        energy_n = openmc.EnergyFilter(x)
        particle_n = openmc.ParticleFilter(['neutron'])
        tally_n = openmc.Tally(name="Neutron Spectrum")
        tally_n.filters = [cell_filter, energy_n, particle_n]
        tally_n.scores = ['flux']
        model.tallies.append(tally_n)
        cwd += 'neutron_spectrum'

    if args.gamma_heating:
        # gamma heating tally 
        settings.survival_biasing = True
        settings.photon_transport = True
        settings.electron_treatment = 'ttb'
        particle_p = openmc.ParticleFilter(['neutron','photon','electron','positron'])
        tally_p = openmc.Tally(name="Gamma heating")
        tally_p.filters = [cell_filter, particle_p]
        tally_p.scores = ['heating']
        model.tallies.append(tally_p)
        cwd = cwd + '-gamma_heating' if cwd else 'gamma_heating'

    if args.reaction_rates:
        # reaction rate tallies 
        tally_1 = openmc.Tally(name="Au-197 (n,gamma)")
        tally_1.nuclides = ['Au197']
        tally_1.filters = [cell_filter]
        tally_1.scores = ['(n,gamma)']

        tally_2 = openmc.Tally(name="In115 (n,n')")
        tally_2.nuclides = ['In115']
        tally_2.filters = [cell_filter]
        tally_2.scores = ['(n,nc)']

        tally_3 = openmc.Tally(name="Al-27 (n,alpha)")
        tally_3.nuclides = ['Al27']
        tally_3.filters = [cell_filter]
        tally_3.scores = ['(n,a)']

        tally_4 = openmc.Tally(name="Nb-93 (n,2n)")
        tally_4.nuclides = ['Nb93']
        tally_4.filters = [cell_filter]
        tally_4.scores = ['(n,2n)']

        tally_5 = openmc.Tally(name="Ni-58 (n,2n)")
        tally_5.nuclides = ['Ni58']
        tally_5.filters = [cell_filter]
        tally_5.scores = ['(n,2n)']

        tally_6 = openmc.Tally(name="Ti-48(n,p)")
        tally_6.nuclides = ['Ti48']
        tally_6.filters = [cell_filter]
        tally_6.scores = ['(n,p)']

        tally_7 = openmc.Tally(name="Ti-49 (n,np)")
        tally_7.nuclides = ['Ti49']
        tally_7.filters = [cell_filter]
        tally_7.scores = ['(n,np)']

        tally_8 = openmc.Tally(name="Ti-47 (n,p)")
        tally_8.nuclides = ['Ti47']
        tally_8.filters = [cell_filter]
        tally_8.scores = ['(n,p)']

        tally_9 = openmc.Tally(name="Ti-48 (n,np)")
        tally_9.nuclides = ['Ti48']
        tally_9.filters = [cell_filter]
        tally_9.scores = ['(n,np)']

        tally_10 = openmc.Tally(name="Fe-56 (n,p)")
        tally_10.nuclides = ['Fe56']
        tally_10.filters = [cell_filter]
        tally_10.scores = ['(n,p)']

        tally_11 = openmc.Tally(name="Zr-90 (n,2n)")
        tally_11.nuclides = ['Zr90']
        tally_11.filters = [cell_filter]
        tally_11.scores = ['(n,2n)']

        tally_12 = openmc.Tally(name="Ni-58 (n,p)")
        tally_12.nuclides = ['Ni58']
        tally_12.filters = [cell_filter]
        tally_12.scores = ['(n,p)']

        tally_13 = openmc.Tally(name="Li-6 (n,tritium)")
        tally_13.nuclides = ['Li6']
        tally_13.filters = [cell_filter]
        tally_13.scores = ['(n,Xt)']

        tally_14 = openmc.Tally(name="Tritium production")
        tally_14.filters = [cell_filter]
        tally_14.scores = ['H3-production']

        tally_15 = openmc.Tally(name="U-235 (n,f)")
        tally_15.nuclides = ['U235']
        tally_15.filters = [cell_filter]
        tally_15.scores = ['fission']

        tally_16 = openmc.Tally(name="Be-9 (n,tritium)")
        tally_16.nuclides = ['Be9']
        tally_16.filters = [cell_filter]
        tally_16.scores = ['(n,Xt)']

        model.tallies.extend([tally_1, tally_2, tally_3, tally_4, tally_5,
                              tally_6, tally_7, tally_8, tally_9, tally_10,
                              tally_11, tally_12, tally_13, tally_14, tally_15,
                              tally_16])

        cwd = cwd + '-rxn_rates' if cwd else 'rxn_rates'

    cwd = usrcwd if usrcwd is not None else cwd
    model.settings = settings

    return model.run(cwd=cwd, threads=args.threads)


if __name__ == "__main__":
    main()
