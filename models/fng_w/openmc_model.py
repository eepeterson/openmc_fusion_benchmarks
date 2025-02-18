#!/usr/bin/env python3
import argparse

import openmc
from openmc_fusion_benchmarks import from_irdff as irdff
from openmc_fusion_benchmarks.neutron_sources import fng_source


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batches", type=int, default=100,
                        help='Number of batches to simulate (int)')
    parser.add_argument("-p", "--particles", type=int,
                        default=int(1e7), help='Number of particles per batch (int)')
    parser.add_argument("-s", "--threads", type=int,
                        help='Number of threads to use in the simulation (int)')
    group = parser.add_argument_group("tallies")
    group.add_argument("-r", "--reaction_rates", action='store_true',
                       default=False, help='Calculate the reaction rates case')
    group.add_argument("-d", "--heating", action='store_true',
                       default=False, help='Calculate the nuclear heating case')

    args = parser.parse_args()

    return args


def main():
    """Analysis of FNG-Tungsten Streaming experiment"""

    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    ############################################################################
    # Define Materials

    mat_1 = openmc.Material(material_id=1)
    mat_1.set_density('g/cm3', 7.954)
    mat_1.add_nuclide('B10', 0.00028, 'wo')
    mat_1.add_nuclide('B11', 0.00287, 'wo')
    mat_1.add_nuclide('C12', 0.039520273803427465, 'wo')
    mat_1.add_nuclide('C13', 0.000479726196572539, 'wo')
    mat_1.add_nuclide('Si28', 0.41, 'wo')
    mat_1.add_nuclide('V50', 0.0003921925373160739, 'wo')
    mat_1.add_nuclide('V51', 0.15960780746268394, 'wo')
    mat_1.add_nuclide('Cr50', 0.724, 'wo')
    mat_1.add_nuclide('Cr52', 14.07, 'wo')
    mat_1.add_nuclide('Cr54', 0.4, 'wo')
    mat_1.add_nuclide('Mn55', 1.14, 'wo')
    mat_1.add_nuclide('Fe54', 4.0, 'wo')
    mat_1.add_nuclide('Fe56', 62.43, 'wo')
    mat_1.add_nuclide('Fe57', 1.48, 'wo')
    mat_1.add_nuclide('Fe58', 0.2, 'wo')
    mat_1.add_nuclide('Co59', 0.14, 'wo')
    mat_1.add_nuclide('Ni58', 7.25, 'wo')
    mat_1.add_nuclide('Ni60', 2.8, 'wo')
    mat_1.add_nuclide('Ni61', 0.13, 'wo')
    mat_1.add_nuclide('Ni62', 0.39, 'wo')
    mat_1.add_nuclide('Ni64', 0.12, 'wo')
    mat_1.add_nuclide('Mo92', 0.2974760972643644, 'wo')
    mat_1.add_nuclide('Mo94', 0.1906159496693566, 'wo')
    mat_1.add_nuclide('Mo95', 0.33284984504744625, 'wo')
    mat_1.add_nuclide('Mo96', 0.3533051306488447, 'wo')
    mat_1.add_nuclide('Mo97', 0.20516502475964524, 'wo')
    mat_1.add_nuclide('Mo98', 0.5254922939545649, 'wo')
    mat_1.add_nuclide('Mo100', 0.21509565865577795, 'wo')
    mat_1.add_nuclide('Cu63', 0.06, 'wo')
    mat_1.add_nuclide('Cu65', 0.03, 'wo')
    #
    mat_2 = openmc.Material(material_id=2)
    mat_2.set_density('g/cm3', 1.0)
    mat_2.add_nuclide('H1', 2.0, 'ao')
    mat_2.add_nuclide('O16', 1.0, 'ao')
    #
    mat_3 = openmc.Material(material_id=3)
    mat_3.set_density('g/cm3', 8.94)
    mat_3.add_nuclide('Cu63', 69.0, 'wo')
    mat_3.add_nuclide('Cu65', 31.0, 'wo')
    #
    mat_4 = openmc.Material(material_id=4)
    mat_4.set_density('atom/b-cm', 4.614e-05)
    mat_4.add_nuclide('N14', 0.788903, 'ao')
    mat_4.add_nuclide('O16', 0.211097, 'ao')
    #
    mat_5 = openmc.Material(material_id=5)
    mat_5.set_density('g/cm3', 18.075)
    mat_5.add_nuclide('W180', 0.0011158467548231742, 'wo')
    mat_5.add_nuclide('W182', 0.24915697002472956, 'wo')
    mat_5.add_nuclide('W183', 0.13528572402024575, 'wo')
    mat_5.add_nuclide('W184', 0.2912528248301395, 'wo')
    mat_5.add_nuclide('W186', 0.2731886343700621, 'wo')
    mat_5.add_nuclide('Ni58', 0.023146, 'wo')
    mat_5.add_nuclide('Ni60', 0.00891582, 'wo')
    mat_5.add_nuclide('Ni61', 0.0003876, 'wo')
    mat_5.add_nuclide('Ni62', 0.00123556, 'wo')
    mat_5.add_nuclide('Ni64', 0.00031484, 'wo')
    mat_5.add_nuclide('Fe54', 0.0009352, 'wo')
    mat_5.add_nuclide('Fe56', 0.0146806, 'wo')
    mat_5.add_nuclide('Fe57', 0.00033904, 'wo')
    mat_5.add_nuclide('Fe58', 4.512e-05, 'wo')
    #
    mat_6 = openmc.Material(material_id=6)
    mat_6.set_density('g/cm3', 0.7)
    mat_6.add_nuclide('Al27', 1.0, 'ao')
    #
    mat_7 = openmc.Material(material_id=7)
    mat_7.set_density('g/cm3', 17.7)
    mat_7.add_nuclide('W180', 0.0010947043952581035, 'wo')
    mat_7.add_nuclide('W182', 0.24443610111899786, 'wo')
    mat_7.add_nuclide('W183', 0.13272241556512532, 'wo')
    mat_7.add_nuclide('W184', 0.2857343502544105, 'wo')
    mat_7.add_nuclide('W186', 0.2680124286662083, 'wo')
    mat_7.add_nuclide('Ni58', 0.0286, 'wo')
    mat_7.add_nuclide('Ni60', 0.011, 'wo')
    mat_7.add_nuclide('Ni61', 0.000479, 'wo')
    mat_7.add_nuclide('Ni62', 0.00153, 'wo')
    mat_7.add_nuclide('Ni64', 0.000389, 'wo')
    mat_7.add_nuclide('Fe54', 0.00152, 'wo')
    mat_7.add_nuclide('Fe56', 0.0239, 'wo')
    mat_7.add_nuclide('Fe57', 0.000551, 'wo')
    mat_7.add_nuclide('Fe58', 7.33e-05, 'wo')
    #
    mat_8 = openmc.Material(material_id=8)
    mat_8.set_density('g/cm3', 2.6)
    mat_8.add_nuclide('H1', 0.005358, 'wo')
    mat_8.add_nuclide('O16', 0.474193, 'wo')
    mat_8.add_nuclide('Na23', 0.016312, 'wo')
    mat_8.add_nuclide('Al27', 0.043491, 'wo')
    mat_8.add_nuclide('Si28', 0.299658, 'wo')
    mat_8.add_nuclide('K39', 0.016757481741937585, 'wo')
    mat_8.add_nuclide('K40', 2.156337426379645e-06, 'wo')
    mat_8.add_nuclide('K41', 0.0012713619206360326, 'wo')
    mat_8.add_nuclide('Fe54', 0.0022582232890827423, 'wo')
    mat_8.add_nuclide('Fe56', 0.036760611511824984, 'wo')
    mat_8.add_nuclide('Fe57', 0.0008641474454698546, 'wo')
    mat_8.add_nuclide('Fe58', 0.00011701775362242353, 'wo')
    # Al27(n,a)Na24
    mat_9 = openmc.Material(material_id=9)
    mat_9.set_density('g/cm3', 2.7)
    mat_9.add_nuclide('Al27', 1.0, 'ao')
    # Fe56(n,p)Mn56
    mat_10 = openmc.Material(material_id=10)
    mat_10.set_density('g/cm3', 7.86)
    mat_10.add_nuclide('Fe54', 0.05845, 'ao')
    mat_10.add_nuclide('Fe56', 0.91754, 'ao')
    mat_10.add_nuclide('Fe57', 0.02119, 'ao')
    mat_10.add_nuclide('Fe58', 0.00282, 'ao')
    # In115(n,n')In115m
    mat_11 = openmc.Material(material_id=11)
    mat_11.set_density('g/cm3', 7.31)
    mat_11.add_nuclide('In113', 0.0429, 'ao')
    mat_11.add_nuclide('In115', 0.9571, 'ao')
    # Nb93(n,2n)Nb92m
    mat_16 = openmc.Material(material_id=16)
    mat_16.set_density('g/cm3', 8.4)
    mat_16.add_nuclide('Nb93', 1.0, 'ao')
    # Ni58(n,2n)Ni57
    mat_26 = openmc.Material(material_id=26)
    mat_26.set_density('g/cm3', 8.9)
    mat_26.add_nuclide('Ni58', 1.0, 'ao')
    # Au197(n,g)Au198
    mat_36 = openmc.Material(material_id=36)
    mat_36.set_density('g/cm3', 19.3)
    mat_36.add_nuclide('Au197', 1.0, 'ao')
    # Zr90(n,2n)Zr89
    mat_46 = openmc.Material(material_id=46)
    mat_46.set_density('g/cm3', 6.5)
    mat_46.add_nuclide('Zr90', 0.5145, 'ao')
    mat_46.add_nuclide('Zr91', 0.1122, 'ao')
    mat_46.add_nuclide('Zr92', 0.1715, 'ao')
    mat_46.add_nuclide('Zr94', 0.1738, 'ao')
    mat_46.add_nuclide('Zr96', 0.028, 'ao')
    # Mn55(n,g)Mn56
    mat_56 = openmc.Material(material_id=56)
    mat_56.set_density('g/cm3', 7.4)
    mat_56.add_nuclide('Mn55', 1.0, 'ao')
    # TLD
    mat_60 = openmc.Material(material_id=60)
    mat_60.set_density('g/cm3', 1.182)
    mat_60.add_nuclide('H1', 0.5333, 'ao')
    mat_60.add_nuclide('C12', 0.3296077026, 'ao')
    mat_60.add_nuclide('C13', 0.0036922973999999995, 'ao')
    mat_60.add_nuclide('O16', 0.1334, 'ao')
    # TLD
    mat_61 = openmc.Material(material_id=61)
    mat_61.set_density('g/cm3', 3.18)
    mat_61.add_nuclide('Ca40', 0.4961920494376943, 'wo')
    mat_61.add_nuclide('Ca42', 0.0034770755162136942, 'wo')
    mat_61.add_nuclide('Ca43', 0.0007428040093963419, 'wo')
    mat_61.add_nuclide('Ca44', 0.011743999254576679, 'wo')
    mat_61.add_nuclide('Ca46', 2.35433917581143e-05, 'wo')
    mat_61.add_nuclide('Ca48', 0.001148528390360927, 'wo')
    mat_61.add_nuclide('F19', 0.486672, 'wo')

    mat_100 = openmc.Material.mix_materials([mat_4, mat_36], [1., 0.0], 'vo')

    # create materials instance
    materials = openmc.Materials([mat_1, mat_2, mat_3, mat_4, mat_5, mat_6, mat_7, mat_8, mat_9,
                                  mat_11, mat_16, mat_26, mat_36, mat_46, mat_56, mat_60, mat_61,
                                  mat_100])

    # %%

    ############################################################################
    # Define Geometry

    # surfaces
    surf_6 = openmc.YPlane(surface_id=6, y0=-1.9)
    surf_1 = openmc.YPlane(surface_id=1, y0=0.0)
    surf_7 = openmc.YCylinder(surface_id=7, x0=0.0, z0=0.0, r=1.5)
    surf_2 = openmc.YPlane(surface_id=2, y0=0.1)
    surf_8 = openmc.YCylinder(surface_id=8, x0=0.0, z0=0.0, r=1.6)
    surf_9 = openmc.YCylinder(surface_id=9, x0=0.0, z0=0.0, r=1.7)
    surf_10 = openmc.YCylinder(surface_id=10, x0=0.0, z0=0.0, r=1.8)
    surf_16 = openmc.Plane(surface_id=16, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=2.0)
    surf_17 = openmc.Plane(surface_id=17, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-2.0)
    surf_30 = openmc.Plane(surface_id=30, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=2.0)
    surf_33 = openmc.Plane(surface_id=33, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-2.0)
    surf_3 = openmc.YPlane(surface_id=3, y0=0.2)
    surf_31 = openmc.Plane(surface_id=31, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=1.9)
    surf_32 = openmc.Plane(surface_id=32, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-1.9)
    surf_4 = openmc.YPlane(surface_id=4, y0=0.35)
    surf_12 = openmc.YPlane(surface_id=12, y0=-13.5)
    surf_13 = openmc.YPlane(surface_id=13, y0=-4.0)
    surf_11 = openmc.YCylinder(surface_id=11, x0=0.0, z0=0.0, r=2.4)
    surf_22 = openmc.Plane(surface_id=22, a=2.6516504294495524,
                           b=-12.450000000000001, c=2.6516504294495524, d=5.01)
    surf_24 = openmc.Plane(surface_id=24, a=2.651650429449554,
                           b=12.450000000000001, c=2.651650429449554, d=8.745)
    surf_42 = openmc.Plane(surface_id=42, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=6.15)
    surf_34 = openmc.Plane(surface_id=34, a=-7.530687219636732,
                           b=8.725608443924892e-16, c=10.076271631908304, d=27.255)
    surf_36 = openmc.Plane(surface_id=36, a=10.076271631908302,
                           b=-6.521244205459656e-16, c=-7.530687219636734, d=27.255)
    surf_23 = openmc.Plane(surface_id=23, a=2.54558441227157, b=-
                           12.450000000000001, c=2.54558441227157, d=2.842499999999999)
    surf_25 = openmc.Plane(surface_id=25, a=2.651650429449554,
                           b=12.450000000000001, c=2.651650429449554, d=7.5)
    surf_35 = openmc.Plane(surface_id=35, a=-7.530687219636731,
                           b=8.725608443924892e-16, c=10.076271631908305, d=28.500000000000004)
    surf_37 = openmc.Plane(surface_id=37, a=10.076271631908304, b=-
                           6.521244205459656e-16, c=-7.530687219636733, d=28.500000000000004)
    surf_26 = openmc.Plane(surface_id=26, a=-2.651650429449554,
                           b=-12.450000000000001, c=-2.651650429449554, d=5.01)
    surf_28 = openmc.Plane(surface_id=28, a=-2.6516504294495524,
                           b=12.450000000000001, c=-2.6516504294495524, d=8.745)
    surf_43 = openmc.Plane(surface_id=43, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-6.15)
    surf_38 = openmc.Plane(surface_id=38, a=-10.076271631908302,
                           b=6.521244205459656e-16, c=7.530687219636734, d=27.255)
    surf_40 = openmc.Plane(surface_id=40, a=7.530687219636732,
                           b=-8.725608443924892e-16, c=-10.076271631908304, d=27.255)
    surf_27 = openmc.Plane(surface_id=27, a=-2.545584412271572, b=-
                           12.450000000000001, c=-2.545584412271572, d=2.842499999999999)
    surf_29 = openmc.Plane(surface_id=29, a=-2.6516504294495524,
                           b=12.450000000000001, c=-2.6516504294495524, d=7.5)
    surf_39 = openmc.Plane(surface_id=39, a=-10.076271631908304,
                           b=6.521244205459656e-16, c=7.530687219636733, d=28.500000000000004)
    surf_41 = openmc.Plane(surface_id=41, a=7.530687219636731, b=-
                           8.725608443924892e-16, c=-10.076271631908305, d=28.500000000000004)
    surf_303 = openmc.ZCylinder(surface_id=303, x0=0.0, y0=-100.0, r=104.9)
    surf_19 = openmc.YCylinder(surface_id=19, x0=0.0, z0=0.0, r=17.7)
    surf_46 = openmc.Plane(surface_id=46, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=10.05)
    surf_47 = openmc.Plane(surface_id=47, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-10.05)
    surf_48 = openmc.Plane(surface_id=48, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=11.45)
    surf_49 = openmc.Plane(surface_id=49, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-11.45)
    surf_57 = openmc.Quadric(surface_id=57, a=0.5000000000000001, b=1.0000000000000002, c=0.5000000000000001, d=-7.850462293418876e-17, e=-
                             7.850462293418876e-17, f=-1.0000000000000002, g=1.1775693440128313e-17, h=-0.30000000000000004, j=1.1775693440128313e-17, k=-3.040000000000001)
    surf_44 = openmc.Plane(surface_id=44, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=7.45)
    surf_15 = openmc.Quadric(surface_id=15, a=0.5000000000000001, b=1.0000000000000002, c=0.5000000000000001, d=-7.850462293418876e-17,
                             e=-7.850462293418876e-17, f=-1.0000000000000002, g=1.1775693440128313e-17, h=-0.30000000000000004, j=1.1775693440128313e-17, k=-1.9375)
    surf_14 = openmc.Quadric(surface_id=14, a=0.5000000000000001, b=1.0000000000000002, c=0.5000000000000001, d=-7.850462293418876e-17,
                             e=-7.850462293418876e-17, f=-1.0000000000000002, g=1.1775693440128313e-17, h=-0.30000000000000004, j=1.1775693440128313e-17, k=-1.2544)
    surf_58 = openmc.Plane(surface_id=58, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=14.35)
    surf_62 = openmc.Plane(surface_id=62, a=6.123233995736766e-17,
                           b=1.0, c=6.123233995736766e-17, d=1.5999999999999999)
    surf_64 = openmc.Plane(
        surface_id=64, a=6.123233995736766e-17, b=1.0, c=6.123233995736766e-17, d=-1.3)
    surf_65 = openmc.Plane(surface_id=65, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=1.45)
    surf_68 = openmc.Plane(surface_id=68, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-1.45)
    surf_60 = openmc.Plane(surface_id=60, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=14.15)
    surf_63 = openmc.Plane(
        surface_id=63, a=6.123233995736766e-17, b=1.0, c=6.123233995736766e-17, d=1.4)
    surf_66 = openmc.Plane(surface_id=66, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=1.25)
    surf_67 = openmc.Plane(surface_id=67, a=-0.7071067811865475,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-1.25)
    surf_50 = openmc.Quadric(surface_id=50, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0,
                             f=-0.0, g=-18.38477631085024, h=0.0, j=-18.38477631085024, k=167.04000000000005)
    surf_51 = openmc.Quadric(surface_id=51, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0,
                             f=-0.0, g=-18.38477631085024, h=0.0, j=-18.38477631085024, k=166.75000000000006)
    surf_20 = openmc.YPlane(surface_id=20, y0=-8.4)
    surf_53 = openmc.Quadric(surface_id=53, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0,
                             f=-0.0, g=-18.38477631085024, h=0.0, j=-18.38477631085024, k=164.59000000000006)
    surf_45 = openmc.Plane(surface_id=45, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-7.45)
    surf_59 = openmc.Plane(surface_id=59, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-14.35)
    surf_61 = openmc.Plane(surface_id=61, a=0.7071067811865476,
                           b=6.123233995736766e-17, c=0.7071067811865476, d=-14.15)
    surf_54 = openmc.Quadric(surface_id=54, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0,
                             f=-0.0, g=18.384776310850242, h=-0.0, j=18.384776310850242, k=167.0400000000001)
    surf_55 = openmc.Quadric(surface_id=55, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0,
                             f=-0.0, g=18.384776310850242, h=-0.0, j=18.384776310850242, k=166.7500000000001)
    surf_56 = openmc.Quadric(surface_id=56, a=1.0, b=0.0, c=1.0, d=-0.0, e=-0.0, f=-
                             0.0, g=18.384776310850242, h=-0.0, j=18.384776310850242, k=164.59000000000012)
    surf_21 = openmc.YPlane(surface_id=21, y0=-14.1, boundary_type="vacuum")
    surf_92 = openmc.XPlane(surface_id=92, x0=-21.0)
    surf_93 = openmc.XPlane(surface_id=93, x0=21.0)
    surf_94 = openmc.ZPlane(surface_id=94, z0=-22.75)
    surf_95 = openmc.ZPlane(surface_id=95, z0=24.1)
    surf_82 = openmc.XPlane(surface_id=82, x0=-23.5)
    surf_98 = openmc.XPlane(surface_id=98, x0=-50.0, boundary_type="vacuum")
    surf_120 = openmc.YPlane(surface_id=120, y0=54.3, boundary_type="vacuum")
    surf_99 = openmc.XPlane(surface_id=99, x0=50.0, boundary_type="vacuum")
    surf_83 = openmc.XPlane(surface_id=83, x0=23.5)
    surf_96 = openmc.ZPlane(surface_id=96, z0=-50.0, boundary_type="vacuum")
    surf_101 = openmc.YPlane(surface_id=101, y0=5.3)
    surf_100 = openmc.ZPlane(surface_id=100, z0=50.0, boundary_type="vacuum")
    surf_84 = openmc.ZPlane(surface_id=84, z0=-15.5)
    surf_85 = openmc.ZPlane(surface_id=85, z0=15.5)
    surf_102 = openmc.YPlane(surface_id=102, y0=7.8)
    surf_201 = openmc.YCylinder(surface_id=201, x0=0.0, z0=0.0, r=0.95)
    surf_202 = openmc.YCylinder(surface_id=202, x0=0.0, z0=0.0, r=2.0)
    surf_203 = openmc.YCylinder(surface_id=203, x0=0.0, z0=0.0, r=3.5)
    surf_204 = openmc.YCylinder(surface_id=204, x0=0.0, z0=0.0, r=5.5)
    surf_86 = openmc.ZPlane(surface_id=86, z0=-3.5)
    surf_87 = openmc.ZPlane(surface_id=87, z0=3.5)
    surf_205 = openmc.YCylinder(surface_id=205, x0=0.0, z0=0.0, r=7.5)
    surf_206 = openmc.YCylinder(surface_id=206, x0=0.0, z0=0.0, r=9.5)
    surf_207 = openmc.YCylinder(surface_id=207, x0=0.0, z0=0.0, r=12.5)
    surf_208 = openmc.YCylinder(surface_id=208, x0=0.0, z0=0.0, r=16.0)
    surf_209 = openmc.YCylinder(surface_id=209, x0=0.0, z0=0.0, r=19.5)
    surf_121 = openmc.YPlane(surface_id=121, y0=10.085)
    surf_500 = openmc.Sphere(surface_id=500, x0=0.0, y0=10.3, z0=0.0, r=0.9)
    surf_124 = openmc.YPlane(surface_id=124, y0=10.515)
    surf_103 = openmc.YPlane(surface_id=103, y0=11.3)
    surf_104 = openmc.YPlane(surface_id=104, y0=13.8)
    surf_105 = openmc.YPlane(surface_id=105, y0=16.5)
    surf_106 = openmc.YPlane(surface_id=106, y0=19.2)
    surf_501 = openmc.Sphere(surface_id=501, x0=0.0, y0=20.3, z0=0.0, r=0.9)
    surf_131 = openmc.YPlane(surface_id=131, y0=20.085)
    surf_134 = openmc.YPlane(surface_id=134, y0=20.515)
    surf_107 = openmc.YPlane(surface_id=107, y0=22.2)
    surf_113 = openmc.YPlane(surface_id=113, y0=25.2)
    surf_114 = openmc.YPlane(surface_id=114, y0=28.5)
    surf_141 = openmc.YPlane(surface_id=141, y0=30.085)
    surf_502 = openmc.Sphere(surface_id=502, x0=0.0, y0=30.3, z0=0.0, r=0.9)
    surf_144 = openmc.YPlane(surface_id=144, y0=30.515)
    surf_115 = openmc.YPlane(surface_id=115, y0=31.8)
    surf_116 = openmc.YPlane(surface_id=116, y0=35.1)
    surf_117 = openmc.YPlane(surface_id=117, y0=38.4)
    surf_151 = openmc.YPlane(surface_id=151, y0=40.085)
    surf_503 = openmc.Sphere(surface_id=503, x0=0.0, y0=40.3, z0=0.0, r=0.9)
    surf_154 = openmc.YPlane(surface_id=154, y0=40.515)
    surf_118 = openmc.YPlane(surface_id=118, y0=42.5)
    surf_119 = openmc.YPlane(surface_id=119, y0=47.0)
    surf_304 = openmc.YCylinder(surface_id=304, x0=0.0, z0=0.0, r=28.0)
    surf_301 = openmc.YCylinder(surface_id=301, x0=0.0, z0=0.0, r=32.5)
    surf_305 = openmc.YCylinder(surface_id=305, x0=0.0, z0=0.0, r=37.0)
    surf_302 = openmc.YCylinder(surface_id=302, x0=0.0, z0=0.0, r=45.0)
    surf_306 = openmc.YCylinder(surface_id=306, x0=0.0, z0=0.0, r=50.0)
    surf_472 = openmc.Plane(surface_id=472, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-810.0)
    surf_473 = openmc.Plane(surface_id=473, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=530.0)
    surf_470 = openmc.Plane(surface_id=470, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-620.0)
    surf_471 = openmc.Plane(surface_id=471, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=620.0)
    surf_474 = openmc.ZPlane(surface_id=474, z0=-456.0)
    surf_475 = openmc.ZPlane(surface_id=475, z0=580.0)
    surf_122 = openmc.YPlane(surface_id=122, y0=10.19)
    surf_125 = openmc.YPlane(surface_id=125, y0=10.29)
    surf_126 = openmc.YPlane(surface_id=126, y0=10.39)
    surf_132 = openmc.YPlane(surface_id=132, y0=20.19)
    surf_135 = openmc.YPlane(surface_id=135, y0=20.29)
    surf_136 = openmc.YPlane(surface_id=136, y0=20.39)
    surf_142 = openmc.YPlane(surface_id=142, y0=30.19)
    surf_145 = openmc.YPlane(surface_id=145, y0=30.29)
    surf_146 = openmc.YPlane(surface_id=146, y0=30.39)
    surf_152 = openmc.YPlane(surface_id=152, y0=40.19)
    surf_155 = openmc.YPlane(surface_id=155, y0=40.29)
    surf_156 = openmc.YPlane(surface_id=156, y0=40.39)
    surf_123 = openmc.YPlane(surface_id=123, y0=10.41)
    surf_133 = openmc.YPlane(surface_id=133, y0=20.41)
    surf_143 = openmc.YPlane(surface_id=143, y0=30.41)
    surf_153 = openmc.YPlane(surface_id=153, y0=40.41)
    surf_402 = openmc.Plane(surface_id=402, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-760.0)
    surf_403 = openmc.Plane(surface_id=403, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=480.0)
    surf_400 = openmc.Plane(surface_id=400, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-570.0)
    surf_401 = openmc.Plane(surface_id=401, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=570.0)
    surf_404 = openmc.ZPlane(surface_id=404, z0=-406.0)
    surf_405 = openmc.ZPlane(surface_id=405, z0=530.0)
    surf_412 = openmc.Plane(surface_id=412, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-762.0)
    surf_413 = openmc.Plane(surface_id=413, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=482.0)
    surf_410 = openmc.Plane(surface_id=410, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-572.0)
    surf_411 = openmc.Plane(surface_id=411, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=572.0)
    surf_414 = openmc.ZPlane(surface_id=414, z0=-408.0)
    surf_415 = openmc.ZPlane(surface_id=415, z0=532.0)
    surf_422 = openmc.Plane(surface_id=422, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-766.0)
    surf_423 = openmc.Plane(surface_id=423, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=486.0)
    surf_420 = openmc.Plane(surface_id=420, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-576.0)
    surf_421 = openmc.Plane(surface_id=421, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=576.0)
    surf_424 = openmc.ZPlane(surface_id=424, z0=-412.0)
    surf_425 = openmc.ZPlane(surface_id=425, z0=536.0)
    surf_432 = openmc.Plane(surface_id=432, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-772.0)
    surf_433 = openmc.Plane(surface_id=433, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=492.0)
    surf_430 = openmc.Plane(surface_id=430, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-582.0)
    surf_431 = openmc.Plane(surface_id=431, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=582.0)
    surf_434 = openmc.ZPlane(surface_id=434, z0=-418.0)
    surf_435 = openmc.ZPlane(surface_id=435, z0=542.0)
    surf_442 = openmc.Plane(surface_id=442, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-780.0)
    surf_443 = openmc.Plane(surface_id=443, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=500.0)
    surf_440 = openmc.Plane(surface_id=440, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-590.0)
    surf_441 = openmc.Plane(surface_id=441, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=590.0)
    surf_444 = openmc.ZPlane(surface_id=444, z0=-426.0)
    surf_445 = openmc.ZPlane(surface_id=445, z0=550.0)
    surf_452 = openmc.Plane(surface_id=452, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-790.0)
    surf_453 = openmc.Plane(surface_id=453, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=510.0)
    surf_450 = openmc.Plane(surface_id=450, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-600.0)
    surf_451 = openmc.Plane(surface_id=451, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=600.0)
    surf_454 = openmc.ZPlane(surface_id=454, z0=-436.0)
    surf_455 = openmc.ZPlane(surface_id=455, z0=560.0)
    surf_462 = openmc.Plane(surface_id=462, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=-800.0)
    surf_463 = openmc.Plane(surface_id=463, a=0.7071067811865476,
                            b=0.7071067811865476, c=6.123233995736766e-17, d=520.0)
    surf_460 = openmc.Plane(surface_id=460, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=-610.0)
    surf_461 = openmc.Plane(surface_id=461, a=0.7071067811865476,
                            b=-0.7071067811865475, c=6.123233995736766e-17, d=610.0)
    surf_464 = openmc.ZPlane(surface_id=464, z0=-446.0)
    surf_465 = openmc.ZPlane(surface_id=465, z0=570.0)

    # regions
    region_1 = (+surf_6 & (+surf_1 | +surf_7) & -surf_2 & -surf_8)
    region_2 = (+surf_6 & +surf_8 & -surf_2 & -surf_9)
    region_3 = (+surf_6 & +surf_9 & -surf_1 & -surf_10)
    region_4 = (+surf_1 & -surf_2 & +surf_9 & -
                surf_16 & +surf_17 & -surf_30 & +surf_33)
    region_5 = (+surf_2 & -surf_3 & -surf_16 & +surf_17 & -surf_31 & +surf_32)
    region_6 = (+surf_2 & -surf_3 & -surf_16 & +surf_17 & -surf_30 & +surf_31)
    region_7 = (+surf_2 & -surf_3 & -surf_16 & +surf_17 & -surf_32 & +surf_33)
    region_8 = (+surf_3 & -surf_4 & -surf_16 & +surf_17 & -surf_30 & +surf_33)
    region_9 = (-surf_6 & +surf_12 & +surf_7 & -surf_10)
    region_10 = (-surf_1 & +surf_12 & -surf_7)
    region_11 = (-surf_6 & +surf_13 & +surf_10 & -surf_11)
    region_12 = (+surf_16 & +surf_22 & +surf_24 & -
                 surf_42 & -surf_34 & -surf_36)
    region_13 = (+surf_16 & +surf_23 & +surf_25 & -surf_42 & -
                 surf_35 & -surf_37 & (-surf_22 | -surf_24 | +surf_34 | +surf_36))
    region_14 = (-surf_17 & +surf_26 & +surf_28 & +
                 surf_43 & -surf_38 & -surf_40)
    region_15 = (-surf_17 & +surf_27 & +surf_29 & +surf_43 & -
                 surf_39 & -surf_41 & (-surf_26 | -surf_28 | +surf_38 | +surf_40))
    region_16 = (+surf_6 & -surf_42 & +surf_43 & -surf_30 & +
                 surf_33 & -surf_25 & -surf_1 & +surf_10 & -surf_29)
    region_17 = (-surf_303 & -surf_42 & +surf_43 & -surf_30 & +
                 surf_33 & -surf_23 & +surf_4 & -surf_27)
    region_18 = (+surf_12 & -surf_13 & +surf_10 & -surf_11)
    region_19 = (+surf_12 & -surf_6 & -surf_19 & +
                 surf_11 & -surf_46 & +surf_47)
    region_20 = (+surf_6 & -surf_303 & -surf_19 & (+surf_42 | -
                 surf_43 | +surf_30 | -surf_33) & -surf_48 & +surf_49 & +surf_57)
    region_21 = (-surf_42 & +surf_16 & -surf_30 & +surf_33 &
                 (+surf_35 | +surf_37) & +surf_23 & +surf_25)
    region_22 = (+surf_42 & -surf_44 & -surf_15 & +surf_14)
    region_23 = (+surf_42 & -surf_44 & -surf_14)
    region_24 = (+surf_42 & -surf_44 & -surf_57 & +surf_15)
    region_25 = (-surf_57 & +surf_44 & -surf_46)
    region_26 = (+surf_46 & -surf_48 & -surf_15 & +surf_14)
    region_27 = (+surf_46 & -surf_48 & -surf_14)
    region_28 = (+surf_46 & -surf_48 & -surf_57 & +surf_15)
    region_29 = (+surf_48 & -surf_58 & -surf_62 & +surf_64 & -
                 surf_65 & +surf_68 & (+surf_60 | +surf_63 | +surf_66 | -surf_67))
    region_30 = (+surf_48 & -surf_60 & -surf_63 & +
                 surf_64 & -surf_66 & +surf_67)
    region_31 = (-surf_64 & -surf_50 & +surf_12)
    region_32 = (-surf_51 & +surf_50 & -surf_64 & +surf_20)
    region_33 = (-surf_53 & +surf_50 & +surf_12 & -surf_20)
    region_34 = (+surf_45 & -surf_43 & -surf_15 & +surf_14)
    region_35 = (+surf_45 & -surf_43 & -surf_14)
    region_36 = (+surf_45 & -surf_43 & -surf_57 & +surf_15)
    region_37 = (-surf_57 & -surf_45 & +surf_47)
    region_38 = (+surf_49 & -surf_47 & -surf_15 & +surf_14)
    region_39 = (+surf_49 & -surf_47 & -surf_14)
    region_40 = (+surf_49 & -surf_47 & -surf_57 & +surf_15)
    region_41 = (-surf_49 & +surf_59 & -surf_62 & +surf_64 & -
                 surf_65 & +surf_68 & (-surf_61 | +surf_63 | +surf_66 | -surf_67))
    region_42 = (-surf_49 & +surf_61 & -surf_63 & +
                 surf_64 & -surf_66 & +surf_67)
    region_43 = (-surf_64 & -surf_54 & +surf_12)
    region_44 = (-surf_55 & +surf_54 & -surf_64 & +surf_20)
    region_45 = (-surf_56 & +surf_54 & +surf_12 & -surf_20)
    region_46 = (+surf_64 & -surf_303 & -surf_19 & +surf_48 &
                 (+surf_58 | +surf_62 | -surf_64 | +surf_65 | -surf_68))
    region_47 = (+surf_64 & -surf_303 & -surf_19 & -surf_49 &
                 (-surf_59 | +surf_62 | -surf_64 | +surf_65 | -surf_68))
    region_48 = (+surf_12 & -surf_64 & -surf_19 & +surf_46 &
                 ((+surf_51 & +surf_20) | (-surf_20 & +surf_53)) & (-surf_6 | +surf_48))
    region_49 = (+surf_12 & -surf_64 & -surf_19 & -surf_47 &
                 ((+surf_55 & +surf_20) | (-surf_20 & +surf_56)) & (-surf_6 | -surf_49))
    region_50 = (-surf_19 & +surf_21 & -surf_12)
    region_51 = (+surf_92 & -surf_93 & -surf_303 & +
                 surf_21 & +surf_94 & -surf_95 & +surf_19)
    region_52 = (-surf_17 & +surf_43 & -surf_30 & +surf_33 &
                 (+surf_39 | +surf_41) & +surf_27 & +surf_29)
    region_53 = (-surf_82 & +surf_98 & -surf_120 & +
                 surf_21 & +surf_94 & -surf_95 & +surf_19)
    region_54 = (-surf_99 & +surf_83 & -surf_120 & +
                 surf_21 & +surf_94 & -surf_95 & +surf_19)
    region_55 = (-surf_94 & +surf_96 & +surf_21 & -
                 surf_101 & -surf_99 & +surf_98)
    region_56 = (-surf_100 & +surf_95 & +surf_21 & -
                 surf_120 & -surf_99 & +surf_98)
    region_57 = (+surf_82 & -surf_92 & +surf_98 & -surf_120 & +
                 surf_101 & +surf_84 & -surf_85 & +surf_19)
    region_58 = (+surf_93 & -surf_83 & -surf_120 & +
                 surf_101 & +surf_84 & -surf_85 & +surf_19)
    region_59 = (+surf_82 & -surf_92 & +surf_98 & -surf_101 & +
                 surf_21 & +surf_84 & -surf_85 & +surf_19)
    region_60 = (+surf_93 & -surf_83 & -surf_101 & +
                 surf_21 & +surf_84 & -surf_85 & +surf_19)
    region_61 = (+surf_82 & -surf_92 & +surf_98 & -surf_101 & +
                 surf_21 & -surf_84 & +surf_94 & +surf_19)
    region_62 = (+surf_82 & -surf_92 & +surf_98 & -surf_101 & +
                 surf_21 & +surf_85 & -surf_95 & +surf_19)
    region_63 = (+surf_93 & -surf_83 & -surf_101 & +
                 surf_21 & +surf_85 & -surf_95 & +surf_19)
    region_64 = (+surf_93 & -surf_83 & +surf_98 & -surf_101 & +
                 surf_21 & -surf_84 & +surf_94 & +surf_19)
    region_101 = (+surf_101 & -surf_102 & -surf_201)
    region_102 = (+surf_101 & -surf_102 & +surf_201 & -surf_202)
    region_103 = (+surf_101 & -surf_102 & +surf_202 & -surf_203)
    region_104 = (+surf_101 & -surf_102 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_105 = (+surf_101 & -surf_102 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_106 = (+surf_101 & -surf_102 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_107 = (+surf_101 & -surf_102 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_108 = (+surf_101 & -surf_102 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_109 = (+surf_101 & -surf_102 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_110 = (+surf_101 & -surf_102 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_111 = (+surf_102 & -surf_121 & -surf_500)
    region_112 = (+surf_124 & -surf_103 & -surf_500)
    region_113 = (+surf_102 & -surf_103 & -surf_201 & +surf_500)
    region_114 = (+surf_102 & -surf_103 & +surf_201 & -surf_202)
    region_115 = (+surf_102 & -surf_103 & +surf_202 & -surf_203)
    region_116 = (+surf_102 & -surf_103 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_117 = (+surf_102 & -surf_103 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_118 = (+surf_102 & -surf_103 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_119 = (+surf_102 & -surf_103 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_120 = (+surf_102 & -surf_103 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_121 = (+surf_102 & -surf_103 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_122 = (+surf_102 & -surf_103 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_124 = (+surf_103 & -surf_104 & -surf_201)
    region_125 = (+surf_103 & -surf_104 & +surf_201 & -surf_202)
    region_126 = (+surf_103 & -surf_104 & +surf_202 & -surf_203)
    region_127 = (+surf_103 & -surf_104 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_128 = (+surf_103 & -surf_104 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_129 = (+surf_103 & -surf_104 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_130 = (+surf_103 & -surf_104 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_131 = (+surf_103 & -surf_104 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_132 = (+surf_103 & -surf_104 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_133 = (+surf_103 & -surf_104 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_136 = (+surf_104 & -surf_105 & -surf_201)
    region_137 = (+surf_104 & -surf_105 & +surf_201 & -surf_202)
    region_138 = (+surf_104 & -surf_105 & +surf_202 & -surf_203)
    region_139 = (+surf_104 & -surf_105 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_140 = (+surf_104 & -surf_105 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_141 = (+surf_104 & -surf_105 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_142 = (+surf_104 & -surf_105 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_143 = (+surf_104 & -surf_105 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_144 = (+surf_104 & -surf_105 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_145 = (+surf_104 & -surf_105 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_147 = (+surf_105 & -surf_106 & -surf_201)
    region_148 = (+surf_105 & -surf_106 & +surf_201 & -surf_202)
    region_149 = (+surf_105 & -surf_106 & +surf_202 & -surf_203)
    region_150 = (+surf_105 & -surf_106 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_151 = (+surf_105 & -surf_106 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_152 = (+surf_105 & -surf_106 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_153 = (+surf_105 & -surf_106 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_154 = (+surf_105 & -surf_106 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_155 = (+surf_105 & -surf_106 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_156 = (+surf_105 & -surf_106 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_157 = (+surf_106 & -surf_501 & -surf_131)
    region_158 = (+surf_134 & -surf_501 & -surf_107)
    region_159 = (+surf_106 & -surf_107 & -surf_201 & +surf_501)
    region_160 = (+surf_106 & -surf_107 & +surf_201 & -surf_202)
    region_161 = (+surf_106 & -surf_107 & +surf_202 & -surf_203)
    region_162 = (+surf_106 & -surf_107 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_163 = (+surf_106 & -surf_107 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_164 = (+surf_106 & -surf_107 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_165 = (+surf_106 & -surf_107 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_166 = (+surf_106 & -surf_107 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_167 = (+surf_106 & -surf_107 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_168 = (+surf_106 & -surf_107 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_228 = (+surf_107 & -surf_113 & -surf_201)
    region_229 = (+surf_107 & -surf_113 & +surf_201 & -surf_202)
    region_230 = (+surf_107 & -surf_113 & +surf_202 & -surf_203)
    region_231 = (+surf_107 & -surf_113 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_232 = (+surf_107 & -surf_113 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_233 = (+surf_107 & -surf_113 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_234 = (+surf_107 & -surf_113 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_235 = (+surf_107 & -surf_113 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_236 = (+surf_107 & -surf_113 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_237 = (+surf_107 & -surf_113 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_240 = (+surf_113 & -surf_114 & -surf_201)
    region_241 = (+surf_113 & -surf_114 & +surf_201 & -surf_202)
    region_242 = (+surf_113 & -surf_114 & +surf_202 & -surf_203)
    region_243 = (+surf_113 & -surf_114 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_244 = (+surf_113 & -surf_114 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_245 = (+surf_113 & -surf_114 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_246 = (+surf_113 & -surf_114 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_247 = (+surf_113 & -surf_114 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_248 = (+surf_113 & -surf_114 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_249 = (+surf_113 & -surf_114 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_250 = (+surf_114 & -surf_141 & -surf_502)
    region_251 = (+surf_144 & -surf_115 & -surf_502)
    region_252 = (+surf_114 & -surf_115 & -surf_201 & +surf_502)
    region_253 = (+surf_114 & -surf_115 & +surf_201 & -surf_202)
    region_254 = (+surf_114 & -surf_115 & +surf_202 & -surf_203)
    region_255 = (+surf_114 & -surf_115 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_256 = (+surf_114 & -surf_115 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_257 = (+surf_114 & -surf_115 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_258 = (+surf_114 & -surf_115 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_259 = (+surf_114 & -surf_115 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_260 = (+surf_114 & -surf_115 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_261 = (+surf_114 & -surf_115 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_263 = (+surf_115 & -surf_116 & -surf_201)
    region_264 = (+surf_115 & -surf_116 & +surf_201 & -surf_202)
    region_265 = (+surf_115 & -surf_116 & +surf_202 & -surf_203)
    region_266 = (+surf_115 & -surf_116 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_267 = (+surf_115 & -surf_116 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_268 = (+surf_115 & -surf_116 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_269 = (+surf_115 & -surf_116 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_270 = (+surf_115 & -surf_116 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_271 = (+surf_115 & -surf_116 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_272 = (+surf_115 & -surf_116 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_274 = (+surf_116 & -surf_117 & -surf_201)
    region_275 = (+surf_116 & -surf_117 & +surf_201 & -surf_202)
    region_276 = (+surf_116 & -surf_117 & +surf_202 & -surf_203)
    region_277 = (+surf_116 & -surf_117 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_278 = (+surf_116 & -surf_117 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_279 = (+surf_116 & -surf_117 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_280 = (+surf_116 & -surf_117 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_281 = (+surf_116 & -surf_117 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_282 = (+surf_116 & -surf_117 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_283 = (+surf_116 & -surf_117 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_284 = (+surf_117 & -surf_151 & -surf_503)
    region_285 = (+surf_154 & -surf_118 & -surf_503)
    region_286 = (+surf_117 & -surf_118 & -surf_201 & +surf_503)
    region_287 = (+surf_117 & -surf_118 & +surf_201 & -surf_202)
    region_288 = (+surf_117 & -surf_118 & +surf_202 & -surf_203)
    region_289 = (+surf_117 & -surf_118 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_290 = (+surf_117 & -surf_118 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_291 = (+surf_117 & -surf_118 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_292 = (+surf_117 & -surf_118 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_293 = (+surf_117 & -surf_118 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_294 = (+surf_117 & -surf_118 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_295 = (+surf_117 & -surf_118 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_297 = (+surf_118 & -surf_119 & -surf_201)
    region_298 = (+surf_118 & -surf_119 & +surf_201 & -surf_202)
    region_299 = (+surf_118 & -surf_119 & +surf_202 & -surf_203)
    region_300 = (+surf_118 & -surf_119 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_301 = (+surf_118 & -surf_119 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_302 = (+surf_118 & -surf_119 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_303 = (+surf_118 & -surf_119 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_304 = (+surf_118 & -surf_119 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_305 = (+surf_118 & -surf_119 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_306 = (+surf_118 & -surf_119 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_309 = (+surf_119 & -surf_120 & -surf_201)
    region_310 = (+surf_119 & -surf_120 & +surf_201 & -surf_202)
    region_311 = (+surf_119 & -surf_120 & +surf_202 & -surf_203)
    region_312 = (+surf_119 & -surf_120 & +surf_203 & -
                  surf_204 & +surf_86 & -surf_87)
    region_313 = (+surf_119 & -surf_120 & +surf_204 & -
                  surf_205 & +surf_86 & -surf_87)
    region_314 = (+surf_119 & -surf_120 & +surf_205 & -
                  surf_206 & +surf_86 & -surf_87)
    region_315 = (+surf_119 & -surf_120 & +surf_206 & -
                  surf_207 & +surf_86 & -surf_87)
    region_316 = (+surf_119 & -surf_120 & +surf_207 & -
                  surf_208 & +surf_86 & -surf_87)
    region_317 = (+surf_119 & -surf_120 & +surf_208 & -
                  surf_209 & +surf_86 & -surf_87)
    region_318 = (+surf_119 & -surf_120 & +surf_209 & +surf_86 & -
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_404 = (-surf_94 & +surf_96 & +surf_101 & -
                  surf_102 & -surf_304 & +surf_98 & -surf_99)
    region_405 = (-surf_94 & +surf_96 & +surf_101 & -surf_102 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_406 = (-surf_94 & +surf_96 & +surf_101 & -surf_102 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_407 = (-surf_94 & +surf_96 & +surf_102 & -
                  surf_103 & -surf_304 & +surf_98 & -surf_99)
    region_408 = (-surf_94 & +surf_96 & +surf_102 & -surf_103 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_409 = (-surf_94 & +surf_96 & +surf_102 & -surf_103 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_411 = (-surf_94 & +surf_96 & +surf_103 & -
                  surf_104 & -surf_304 & +surf_98 & -surf_99)
    region_412 = (-surf_94 & +surf_96 & +surf_103 & -surf_104 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_413 = (-surf_94 & +surf_96 & +surf_103 & -surf_104 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_414 = (-surf_94 & +surf_96 & +surf_104 & -
                  surf_105 & -surf_304 & +surf_98 & -surf_99)
    region_415 = (-surf_94 & +surf_96 & +surf_104 & -surf_105 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_416 = (-surf_94 & +surf_96 & +surf_104 & -surf_105 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_417 = (-surf_94 & +surf_96 & +surf_105 & -
                  surf_106 & -surf_304 & +surf_98 & -surf_99)
    region_418 = (-surf_94 & +surf_96 & +surf_105 & -surf_106 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_419 = (-surf_94 & +surf_96 & +surf_105 & -surf_106 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_421 = (-surf_94 & +surf_96 & +surf_106 & -
                  surf_107 & -surf_304 & +surf_98 & -surf_99)
    region_422 = (-surf_94 & +surf_96 & +surf_106 & -surf_107 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_423 = (-surf_94 & +surf_96 & +surf_106 & -surf_107 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_424 = (-surf_94 & +surf_96 & +surf_107 & -
                  surf_113 & -surf_304 & +surf_98 & -surf_99)
    region_425 = (-surf_94 & +surf_96 & +surf_107 & -surf_113 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_426 = (-surf_94 & +surf_96 & +surf_107 & -surf_113 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_427 = (-surf_94 & +surf_96 & +surf_113 & -
                  surf_114 & -surf_304 & +surf_98 & -surf_99)
    region_428 = (-surf_94 & +surf_96 & +surf_113 & -surf_114 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_429 = (-surf_94 & +surf_96 & +surf_113 & -surf_114 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_431 = (-surf_94 & +surf_96 & +surf_114 & -
                  surf_115 & -surf_304 & +surf_98 & -surf_99)
    region_432 = (-surf_94 & +surf_96 & +surf_114 & -surf_115 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_433 = (-surf_94 & +surf_96 & +surf_114 & -surf_115 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_434 = (-surf_94 & +surf_96 & +surf_115 & -
                  surf_116 & -surf_304 & +surf_98 & -surf_99)
    region_435 = (-surf_94 & +surf_96 & +surf_115 & -surf_116 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_436 = (-surf_94 & +surf_96 & +surf_115 & -surf_116 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_437 = (-surf_94 & +surf_96 & +surf_116 & -
                  surf_117 & -surf_304 & +surf_98 & -surf_99)
    region_438 = (-surf_94 & +surf_96 & +surf_116 & -surf_117 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_439 = (-surf_94 & +surf_96 & +surf_116 & -surf_117 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_441 = (-surf_94 & +surf_96 & +surf_117 & -
                  surf_118 & -surf_304 & +surf_98 & -surf_99)
    region_442 = (-surf_94 & +surf_96 & +surf_117 & -surf_118 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_443 = (-surf_94 & +surf_96 & +surf_117 & -surf_118 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_444 = (-surf_94 & +surf_96 & +surf_118 & -
                  surf_119 & -surf_304 & +surf_98 & -surf_99)
    region_445 = (-surf_94 & +surf_96 & +surf_118 & -surf_119 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_446 = (-surf_94 & +surf_96 & +surf_118 & -surf_119 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_447 = (-surf_94 & +surf_96 & +surf_119 & -
                  surf_120 & -surf_304 & +surf_98 & -surf_99)
    region_448 = (-surf_94 & +surf_96 & +surf_119 & -surf_120 & +
                  surf_304 & -surf_301 & +surf_98 & -surf_99)
    region_449 = (-surf_94 & +surf_96 & +surf_119 & -surf_120 & +
                  surf_301 & -surf_305 & +surf_98 & -surf_99)
    region_504 = (-surf_94 & +surf_96 & +surf_101 & -surf_102 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_505 = (-surf_94 & +surf_96 & +surf_101 & -surf_102 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_506 = (-surf_94 & +surf_96 & +surf_101 & -
                  surf_102 & +surf_306 & +surf_98 & -surf_99)
    region_507 = (-surf_94 & +surf_96 & +surf_102 & -surf_103 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_508 = (-surf_94 & +surf_96 & +surf_102 & -surf_103 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_509 = (-surf_94 & +surf_96 & +surf_102 & -
                  surf_103 & +surf_306 & +surf_98 & -surf_99)
    region_511 = (-surf_94 & +surf_96 & +surf_103 & -surf_104 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_512 = (-surf_94 & +surf_96 & +surf_103 & -surf_104 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_513 = (-surf_94 & +surf_96 & +surf_103 & -
                  surf_104 & +surf_306 & +surf_98 & -surf_99)
    region_514 = (-surf_94 & +surf_96 & +surf_104 & -surf_105 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_515 = (-surf_94 & +surf_96 & +surf_104 & -surf_105 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_516 = (-surf_94 & +surf_96 & +surf_104 & -
                  surf_105 & +surf_306 & +surf_98 & -surf_99)
    region_517 = (-surf_94 & +surf_96 & +surf_105 & -surf_106 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_518 = (-surf_94 & +surf_96 & +surf_105 & -surf_106 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_519 = (-surf_94 & +surf_96 & +surf_105 & -
                  surf_106 & +surf_306 & +surf_98 & -surf_99)
    region_521 = (-surf_94 & +surf_96 & +surf_106 & -surf_107 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_522 = (-surf_94 & +surf_96 & +surf_106 & -surf_107 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_523 = (-surf_94 & +surf_96 & +surf_106 & -
                  surf_107 & +surf_306 & +surf_98 & -surf_99)
    region_524 = (-surf_94 & +surf_96 & +surf_107 & -surf_113 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_525 = (-surf_94 & +surf_96 & +surf_107 & -surf_113 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_526 = (-surf_94 & +surf_96 & +surf_107 & -
                  surf_113 & +surf_306 & +surf_98 & -surf_99)
    region_527 = (-surf_94 & +surf_96 & +surf_113 & -surf_114 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_528 = (-surf_94 & +surf_96 & +surf_113 & -surf_114 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_529 = (-surf_94 & +surf_96 & +surf_113 & -
                  surf_114 & +surf_306 & +surf_98 & -surf_99)
    region_531 = (-surf_94 & +surf_96 & +surf_114 & -surf_115 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_532 = (-surf_94 & +surf_96 & +surf_114 & -surf_115 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_533 = (-surf_94 & +surf_96 & +surf_114 & -
                  surf_115 & +surf_306 & +surf_98 & -surf_99)
    region_534 = (-surf_94 & +surf_96 & +surf_115 & -surf_116 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_535 = (-surf_94 & +surf_96 & +surf_115 & -surf_116 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_536 = (-surf_94 & +surf_96 & +surf_115 & -
                  surf_116 & +surf_306 & +surf_98 & -surf_99)
    region_537 = (-surf_94 & +surf_96 & +surf_116 & -surf_117 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_538 = (-surf_94 & +surf_96 & +surf_116 & -surf_117 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_539 = (-surf_94 & +surf_96 & +surf_116 & -
                  surf_117 & +surf_306 & +surf_98 & -surf_99)
    region_541 = (-surf_94 & +surf_96 & +surf_117 & -surf_118 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_542 = (-surf_94 & +surf_96 & +surf_117 & -surf_118 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_543 = (-surf_94 & +surf_96 & +surf_117 & -
                  surf_118 & +surf_306 & +surf_98 & -surf_99)
    region_544 = (-surf_94 & +surf_96 & +surf_118 & -surf_119 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_545 = (-surf_94 & +surf_96 & +surf_118 & -surf_119 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_546 = (-surf_94 & +surf_96 & +surf_118 & -
                  surf_119 & +surf_306 & +surf_98 & -surf_99)
    region_547 = (-surf_94 & +surf_96 & +surf_119 & -surf_120 & +
                  surf_305 & -surf_302 & +surf_98 & -surf_99)
    region_548 = (-surf_94 & +surf_96 & +surf_119 & -surf_120 & +
                  surf_302 & -surf_306 & +surf_98 & -surf_99)
    region_549 = (-surf_94 & +surf_96 & +surf_119 & -
                  surf_120 & +surf_306 & +surf_98 & -surf_99)
    region_600 = (-surf_472 | +surf_473 | -surf_470 | +
                  surf_471 | -surf_474 | +surf_475)
    region_601 = (+surf_121 & -surf_122 & -surf_500)
    region_602 = (+surf_122 & -surf_125 & -surf_500)
    region_611 = (+surf_125 & -surf_126 & -surf_500)
    region_612 = (+surf_131 & -surf_132 & -surf_501)
    region_613 = (+surf_132 & -surf_135 & -surf_501)
    region_614 = (+surf_135 & -surf_136 & -surf_501)
    region_621 = (+surf_141 & -surf_142 & -surf_502)
    region_622 = (+surf_142 & -surf_145 & -surf_502)
    region_623 = (+surf_145 & -surf_146 & -surf_502)
    region_631 = (+surf_151 & -surf_152 & -surf_503)
    region_632 = (+surf_152 & -surf_155 & -surf_503)
    region_633 = (+surf_155 & -surf_156 & -surf_503)
    region_641 = (+surf_126 & -surf_123 & -surf_500)
    region_642 = (+surf_123 & -surf_124 & -surf_500)
    region_643 = (+surf_136 & -surf_133 & -surf_501)
    region_644 = (+surf_133 & -surf_134 & -surf_501)
    region_645 = (+surf_146 & -surf_143 & -surf_502)
    region_646 = (+surf_143 & -surf_144 & -surf_502)
    region_647 = (+surf_156 & -surf_153 & -surf_503)
    region_648 = (+surf_153 & -surf_154 & -surf_503)
    region_654 = (+surf_402 & -surf_403 & +surf_400 & -surf_401 & +surf_404 & -
                  surf_405 & (+surf_120 | -surf_98 | +surf_99 | -surf_96 | +surf_100 | -surf_21))
    region_655 = ((-surf_402 | +surf_403 | -surf_400 | +surf_401 | -surf_404 | +surf_405)
                  & +surf_412 & -surf_413 & +surf_410 & -surf_411 & +surf_414 & -surf_415)
    region_656 = ((-surf_412 | +surf_413 | -surf_410 | +surf_411 | -surf_414 | +surf_415)
                  & +surf_422 & -surf_423 & +surf_420 & -surf_421 & +surf_424 & -surf_425)
    region_657 = ((-surf_422 | +surf_423 | -surf_420 | +surf_421 | -surf_424 | +surf_425)
                  & +surf_432 & -surf_433 & +surf_430 & -surf_431 & +surf_434 & -surf_435)
    region_658 = ((-surf_432 | +surf_433 | -surf_430 | +surf_431 | -surf_434 | +surf_435)
                  & +surf_442 & -surf_443 & +surf_440 & -surf_441 & +surf_444 & -surf_445)
    region_659 = ((-surf_442 | +surf_443 | -surf_440 | +surf_441 | -surf_444 | +surf_445)
                  & +surf_452 & -surf_453 & +surf_450 & -surf_451 & +surf_454 & -surf_455)
    region_660 = ((-surf_452 | +surf_453 | -surf_450 | +surf_451 | -surf_454 | +surf_455)
                  & +surf_462 & -surf_463 & +surf_460 & -surf_461 & +surf_464 & -surf_465)
    region_661 = ((-surf_462 | +surf_463 | -surf_460 | +surf_461 | -surf_464 | +surf_465)
                  & +surf_472 & -surf_473 & +surf_470 & -surf_471 & +surf_474 & -surf_475)
    region_701 = (+surf_101 & -surf_102 & +surf_203 & -surf_204 & +surf_87)
    region_702 = (+surf_101 & -surf_102 & +surf_204 & -surf_205 & +surf_87)
    region_703 = (+surf_101 & -surf_102 & +surf_205 & -surf_206 & +surf_87)
    region_704 = (+surf_101 & -surf_102 & +surf_206 & -surf_207 & +surf_87)
    region_705 = (+surf_101 & -surf_102 & +surf_207 & -surf_208 & +surf_87)
    region_706 = (+surf_101 & -surf_102 & +surf_208 & -surf_209 & +surf_87)
    region_707 = (+surf_101 & -surf_102 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_708 = (+surf_102 & -surf_103 & +surf_203 & -surf_204 & +surf_87)
    region_709 = (+surf_102 & -surf_103 & +surf_204 & -surf_205 & +surf_87)
    region_710 = (+surf_102 & -surf_103 & +surf_205 & -surf_206 & +surf_87)
    region_711 = (+surf_102 & -surf_103 & +surf_206 & -surf_207 & +surf_87)
    region_712 = (+surf_102 & -surf_103 & +surf_207 & -surf_208 & +surf_87)
    region_713 = (+surf_102 & -surf_103 & +surf_208 & -surf_209 & +surf_87)
    region_714 = (+surf_102 & -surf_103 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_715 = (+surf_103 & -surf_104 & +surf_203 & -surf_204 & +surf_87)
    region_716 = (+surf_103 & -surf_104 & +surf_204 & -surf_205 & +surf_87)
    region_717 = (+surf_103 & -surf_104 & +surf_205 & -surf_206 & +surf_87)
    region_718 = (+surf_103 & -surf_104 & +surf_206 & -surf_207 & +surf_87)
    region_719 = (+surf_103 & -surf_104 & +surf_207 & -surf_208 & +surf_87)
    region_720 = (+surf_103 & -surf_104 & +surf_208 & -surf_209 & +surf_87)
    region_721 = (+surf_103 & -surf_104 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_722 = (+surf_104 & -surf_105 & +surf_203 & -surf_204 & +surf_87)
    region_723 = (+surf_104 & -surf_105 & +surf_204 & -surf_205 & +surf_87)
    region_724 = (+surf_104 & -surf_105 & +surf_205 & -surf_206 & +surf_87)
    region_725 = (+surf_104 & -surf_105 & +surf_206 & -surf_207 & +surf_87)
    region_726 = (+surf_104 & -surf_105 & +surf_207 & -surf_208 & +surf_87)
    region_727 = (+surf_104 & -surf_105 & +surf_208 & -surf_209 & +surf_87)
    region_728 = (+surf_104 & -surf_105 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_729 = (+surf_105 & -surf_106 & +surf_203 & -surf_204 & +surf_87)
    region_730 = (+surf_105 & -surf_106 & +surf_204 & -surf_205 & +surf_87)
    region_731 = (+surf_105 & -surf_106 & +surf_205 & -surf_206 & +surf_87)
    region_732 = (+surf_105 & -surf_106 & +surf_206 & -surf_207 & +surf_87)
    region_733 = (+surf_105 & -surf_106 & +surf_207 & -surf_208 & +surf_87)
    region_734 = (+surf_105 & -surf_106 & +surf_208 & -surf_209 & +surf_87)
    region_735 = (+surf_105 & -surf_106 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_736 = (+surf_106 & -surf_107 & +surf_203 & -surf_204 & +surf_87)
    region_737 = (+surf_106 & -surf_107 & +surf_204 & -surf_205 & +surf_87)
    region_738 = (+surf_106 & -surf_107 & +surf_205 & -surf_206 & +surf_87)
    region_739 = (+surf_106 & -surf_107 & +surf_206 & -surf_207 & +surf_87)
    region_740 = (+surf_106 & -surf_107 & +surf_207 & -surf_208 & +surf_87)
    region_741 = (+surf_106 & -surf_107 & +surf_208 & -surf_209 & +surf_87)
    region_742 = (+surf_106 & -surf_107 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_743 = (+surf_107 & -surf_113 & +surf_203 & -surf_204 & +surf_87)
    region_744 = (+surf_107 & -surf_113 & +surf_204 & -surf_205 & +surf_87)
    region_745 = (+surf_107 & -surf_113 & +surf_205 & -surf_206 & +surf_87)
    region_746 = (+surf_107 & -surf_113 & +surf_206 & -surf_207 & +surf_87)
    region_747 = (+surf_107 & -surf_113 & +surf_207 & -surf_208 & +surf_87)
    region_748 = (+surf_107 & -surf_113 & +surf_208 & -surf_209 & +surf_87)
    region_749 = (+surf_107 & -surf_113 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_750 = (+surf_113 & -surf_114 & +surf_203 & -surf_204 & +surf_87)
    region_751 = (+surf_113 & -surf_114 & +surf_204 & -surf_205 & +surf_87)
    region_752 = (+surf_113 & -surf_114 & +surf_205 & -surf_206 & +surf_87)
    region_753 = (+surf_113 & -surf_114 & +surf_206 & -surf_207 & +surf_87)
    region_754 = (+surf_113 & -surf_114 & +surf_207 & -surf_208 & +surf_87)
    region_755 = (+surf_113 & -surf_114 & +surf_208 & -surf_209 & +surf_87)
    region_756 = (+surf_113 & -surf_114 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_757 = (+surf_114 & -surf_115 & +surf_203 & -surf_204 & +surf_87)
    region_758 = (+surf_114 & -surf_115 & +surf_204 & -surf_205 & +surf_87)
    region_759 = (+surf_114 & -surf_115 & +surf_205 & -surf_206 & +surf_87)
    region_760 = (+surf_114 & -surf_115 & +surf_206 & -surf_207 & +surf_87)
    region_761 = (+surf_114 & -surf_115 & +surf_207 & -surf_208 & +surf_87)
    region_762 = (+surf_114 & -surf_115 & +surf_208 & -surf_209 & +surf_87)
    region_763 = (+surf_114 & -surf_115 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_764 = (+surf_115 & -surf_116 & +surf_203 & -surf_204 & +surf_87)
    region_765 = (+surf_115 & -surf_116 & +surf_204 & -surf_205 & +surf_87)
    region_766 = (+surf_115 & -surf_116 & +surf_205 & -surf_206 & +surf_87)
    region_767 = (+surf_115 & -surf_116 & +surf_206 & -surf_207 & +surf_87)
    region_768 = (+surf_115 & -surf_116 & +surf_207 & -surf_208 & +surf_87)
    region_769 = (+surf_115 & -surf_116 & +surf_208 & -surf_209 & +surf_87)
    region_770 = (+surf_115 & -surf_116 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_771 = (+surf_116 & -surf_117 & +surf_203 & -surf_204 & +surf_87)
    region_772 = (+surf_116 & -surf_117 & +surf_204 & -surf_205 & +surf_87)
    region_773 = (+surf_116 & -surf_117 & +surf_205 & -surf_206 & +surf_87)
    region_774 = (+surf_116 & -surf_117 & +surf_206 & -surf_207 & +surf_87)
    region_775 = (+surf_116 & -surf_117 & +surf_207 & -surf_208 & +surf_87)
    region_776 = (+surf_116 & -surf_117 & +surf_208 & -surf_209 & +surf_87)
    region_777 = (+surf_116 & -surf_117 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_778 = (+surf_117 & -surf_118 & +surf_203 & -surf_204 & +surf_87)
    region_779 = (+surf_117 & -surf_118 & +surf_204 & -surf_205 & +surf_87)
    region_780 = (+surf_117 & -surf_118 & +surf_205 & -surf_206 & +surf_87)
    region_781 = (+surf_117 & -surf_118 & +surf_206 & -surf_207 & +surf_87)
    region_782 = (+surf_117 & -surf_118 & +surf_207 & -surf_208 & +surf_87)
    region_783 = (+surf_117 & -surf_118 & +surf_208 & -surf_209 & +surf_87)
    region_784 = (+surf_117 & -surf_118 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_785 = (+surf_118 & -surf_119 & +surf_203 & -surf_204 & +surf_87)
    region_786 = (+surf_118 & -surf_119 & +surf_204 & -surf_205 & +surf_87)
    region_787 = (+surf_118 & -surf_119 & +surf_205 & -surf_206 & +surf_87)
    region_788 = (+surf_118 & -surf_119 & +surf_206 & -surf_207 & +surf_87)
    region_789 = (+surf_118 & -surf_119 & +surf_207 & -surf_208 & +surf_87)
    region_790 = (+surf_118 & -surf_119 & +surf_208 & -surf_209 & +surf_87)
    region_791 = (+surf_118 & -surf_119 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_792 = (+surf_119 & -surf_120 & +surf_203 & -surf_204 & +surf_87)
    region_793 = (+surf_119 & -surf_120 & +surf_204 & -surf_205 & +surf_87)
    region_794 = (+surf_119 & -surf_120 & +surf_205 & -surf_206 & +surf_87)
    region_795 = (+surf_119 & -surf_120 & +surf_206 & -surf_207 & +surf_87)
    region_796 = (+surf_119 & -surf_120 & +surf_207 & -surf_208 & +surf_87)
    region_797 = (+surf_119 & -surf_120 & +surf_208 & -surf_209 & +surf_87)
    region_798 = (+surf_119 & -surf_120 & +surf_209 & +
                  surf_87 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_801 = (+surf_101 & -surf_102 & +surf_203 & -surf_204 & -surf_86)
    region_802 = (+surf_101 & -surf_102 & +surf_204 & -surf_205 & -surf_86)
    region_803 = (+surf_101 & -surf_102 & +surf_205 & -surf_206 & -surf_86)
    region_804 = (+surf_101 & -surf_102 & +surf_206 & -surf_207 & -surf_86)
    region_805 = (+surf_101 & -surf_102 & +surf_207 & -surf_208 & -surf_86)
    region_806 = (+surf_101 & -surf_102 & +surf_208 & -surf_209 & -surf_86)
    region_807 = (+surf_101 & -surf_102 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_808 = (+surf_102 & -surf_103 & +surf_203 & -surf_204 & -surf_86)
    region_809 = (+surf_102 & -surf_103 & +surf_204 & -surf_205 & -surf_86)
    region_810 = (+surf_102 & -surf_103 & +surf_205 & -surf_206 & -surf_86)
    region_811 = (+surf_102 & -surf_103 & +surf_206 & -surf_207 & -surf_86)
    region_812 = (+surf_102 & -surf_103 & +surf_207 & -surf_208 & -surf_86)
    region_813 = (+surf_102 & -surf_103 & +surf_208 & -surf_209 & -surf_86)
    region_814 = (+surf_102 & -surf_103 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_815 = (+surf_103 & -surf_104 & +surf_203 & -surf_204 & -surf_86)
    region_816 = (+surf_103 & -surf_104 & +surf_204 & -surf_205 & -surf_86)
    region_817 = (+surf_103 & -surf_104 & +surf_205 & -surf_206 & -surf_86)
    region_818 = (+surf_103 & -surf_104 & +surf_206 & -surf_207 & -surf_86)
    region_819 = (+surf_103 & -surf_104 & +surf_207 & -surf_208 & -surf_86)
    region_820 = (+surf_103 & -surf_104 & +surf_208 & -surf_209 & -surf_86)
    region_821 = (+surf_103 & -surf_104 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_822 = (+surf_104 & -surf_105 & +surf_203 & -surf_204 & -surf_86)
    region_823 = (+surf_104 & -surf_105 & +surf_204 & -surf_205 & -surf_86)
    region_824 = (+surf_104 & -surf_105 & +surf_205 & -surf_206 & -surf_86)
    region_825 = (+surf_104 & -surf_105 & +surf_206 & -surf_207 & -surf_86)
    region_826 = (+surf_104 & -surf_105 & +surf_207 & -surf_208 & -surf_86)
    region_827 = (+surf_104 & -surf_105 & +surf_208 & -surf_209 & -surf_86)
    region_828 = (+surf_104 & -surf_105 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_829 = (+surf_105 & -surf_106 & +surf_203 & -surf_204 & -surf_86)
    region_830 = (+surf_105 & -surf_106 & +surf_204 & -surf_205 & -surf_86)
    region_831 = (+surf_105 & -surf_106 & +surf_205 & -surf_206 & -surf_86)
    region_832 = (+surf_105 & -surf_106 & +surf_206 & -surf_207 & -surf_86)
    region_833 = (+surf_105 & -surf_106 & +surf_207 & -surf_208 & -surf_86)
    region_834 = (+surf_105 & -surf_106 & +surf_208 & -surf_209 & -surf_86)
    region_835 = (+surf_105 & -surf_106 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_836 = (+surf_106 & -surf_107 & +surf_203 & -surf_204 & -surf_86)
    region_837 = (+surf_106 & -surf_107 & +surf_204 & -surf_205 & -surf_86)
    region_838 = (+surf_106 & -surf_107 & +surf_205 & -surf_206 & -surf_86)
    region_839 = (+surf_106 & -surf_107 & +surf_206 & -surf_207 & -surf_86)
    region_840 = (+surf_106 & -surf_107 & +surf_207 & -surf_208 & -surf_86)
    region_841 = (+surf_106 & -surf_107 & +surf_208 & -surf_209 & -surf_86)
    region_842 = (+surf_106 & -surf_107 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_843 = (+surf_107 & -surf_113 & +surf_203 & -surf_204 & -surf_86)
    region_844 = (+surf_107 & -surf_113 & +surf_204 & -surf_205 & -surf_86)
    region_845 = (+surf_107 & -surf_113 & +surf_205 & -surf_206 & -surf_86)
    region_846 = (+surf_107 & -surf_113 & +surf_206 & -surf_207 & -surf_86)
    region_847 = (+surf_107 & -surf_113 & +surf_207 & -surf_208 & -surf_86)
    region_848 = (+surf_107 & -surf_113 & +surf_208 & -surf_209 & -surf_86)
    region_849 = (+surf_107 & -surf_113 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_850 = (+surf_113 & -surf_114 & +surf_203 & -surf_204 & -surf_86)
    region_851 = (+surf_113 & -surf_114 & +surf_204 & -surf_205 & -surf_86)
    region_852 = (+surf_113 & -surf_114 & +surf_205 & -surf_206 & -surf_86)
    region_853 = (+surf_113 & -surf_114 & +surf_206 & -surf_207 & -surf_86)
    region_854 = (+surf_113 & -surf_114 & +surf_207 & -surf_208 & -surf_86)
    region_855 = (+surf_113 & -surf_114 & +surf_208 & -surf_209 & -surf_86)
    region_856 = (+surf_113 & -surf_114 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_857 = (+surf_114 & -surf_115 & +surf_203 & -surf_204 & -surf_86)
    region_858 = (+surf_114 & -surf_115 & +surf_204 & -surf_205 & -surf_86)
    region_859 = (+surf_114 & -surf_115 & +surf_205 & -surf_206 & -surf_86)
    region_860 = (+surf_114 & -surf_115 & +surf_206 & -surf_207 & -surf_86)
    region_861 = (+surf_114 & -surf_115 & +surf_207 & -surf_208 & -surf_86)
    region_862 = (+surf_114 & -surf_115 & +surf_208 & -surf_209 & -surf_86)
    region_863 = (+surf_114 & -surf_115 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_864 = (+surf_115 & -surf_116 & +surf_203 & -surf_204 & -surf_86)
    region_865 = (+surf_115 & -surf_116 & +surf_204 & -surf_205 & -surf_86)
    region_866 = (+surf_115 & -surf_116 & +surf_205 & -surf_206 & -surf_86)
    region_867 = (+surf_115 & -surf_116 & +surf_206 & -surf_207 & -surf_86)
    region_868 = (+surf_115 & -surf_116 & +surf_207 & -surf_208 & -surf_86)
    region_869 = (+surf_115 & -surf_116 & +surf_208 & -surf_209 & -surf_86)
    region_870 = (+surf_115 & -surf_116 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_871 = (+surf_116 & -surf_117 & +surf_203 & -surf_204 & -surf_86)
    region_872 = (+surf_116 & -surf_117 & +surf_204 & -surf_205 & -surf_86)
    region_873 = (+surf_116 & -surf_117 & +surf_205 & -surf_206 & -surf_86)
    region_874 = (+surf_116 & -surf_117 & +surf_206 & -surf_207 & -surf_86)
    region_875 = (+surf_116 & -surf_117 & +surf_207 & -surf_208 & -surf_86)
    region_876 = (+surf_116 & -surf_117 & +surf_208 & -surf_209 & -surf_86)
    region_877 = (+surf_116 & -surf_117 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_878 = (+surf_117 & -surf_118 & +surf_203 & -surf_204 & -surf_86)
    region_879 = (+surf_117 & -surf_118 & +surf_204 & -surf_205 & -surf_86)
    region_880 = (+surf_117 & -surf_118 & +surf_205 & -surf_206 & -surf_86)
    region_881 = (+surf_117 & -surf_118 & +surf_206 & -surf_207 & -surf_86)
    region_882 = (+surf_117 & -surf_118 & +surf_207 & -surf_208 & -surf_86)
    region_883 = (+surf_117 & -surf_118 & +surf_208 & -surf_209 & -surf_86)
    region_884 = (+surf_117 & -surf_118 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_885 = (+surf_118 & -surf_119 & +surf_203 & -surf_204 & -surf_86)
    region_886 = (+surf_118 & -surf_119 & +surf_204 & -surf_205 & -surf_86)
    region_887 = (+surf_118 & -surf_119 & +surf_205 & -surf_206 & -surf_86)
    region_888 = (+surf_118 & -surf_119 & +surf_206 & -surf_207 & -surf_86)
    region_889 = (+surf_118 & -surf_119 & +surf_207 & -surf_208 & -surf_86)
    region_890 = (+surf_118 & -surf_119 & +surf_208 & -surf_209 & -surf_86)
    region_891 = (+surf_118 & -surf_119 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_892 = (+surf_119 & -surf_120 & +surf_203 & -surf_204 & -surf_86)
    region_893 = (+surf_119 & -surf_120 & +surf_204 & -surf_205 & -surf_86)
    region_894 = (+surf_119 & -surf_120 & +surf_205 & -surf_206 & -surf_86)
    region_895 = (+surf_119 & -surf_120 & +surf_206 & -surf_207 & -surf_86)
    region_896 = (+surf_119 & -surf_120 & +surf_207 & -surf_208 & -surf_86)
    region_897 = (+surf_119 & -surf_120 & +surf_208 & -surf_209 & -surf_86)
    region_898 = (+surf_119 & -surf_120 & +surf_209 & -
                  surf_86 & +surf_92 & -surf_93 & +surf_94 & -surf_95)
    region_900 = (+surf_101 & -surf_102 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_901 = (+surf_101 & -surf_102 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_902 = (+surf_101 & -surf_102 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_903 = (+surf_101 & -surf_102 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_904 = (+surf_102 & -surf_103 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_905 = (+surf_102 & -surf_103 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_906 = (+surf_102 & -surf_103 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_907 = (+surf_102 & -surf_103 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_908 = (+surf_103 & -surf_104 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_909 = (+surf_103 & -surf_104 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_910 = (+surf_103 & -surf_104 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_911 = (+surf_103 & -surf_104 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_912 = (+surf_104 & -surf_105 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_913 = (+surf_104 & -surf_105 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_914 = (+surf_104 & -surf_105 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_915 = (+surf_104 & -surf_105 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_916 = (+surf_105 & -surf_106 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_917 = (+surf_105 & -surf_106 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_918 = (+surf_105 & -surf_106 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_919 = (+surf_105 & -surf_106 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_920 = (+surf_106 & -surf_107 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_921 = (+surf_106 & -surf_107 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_922 = (+surf_106 & -surf_107 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_923 = (+surf_106 & -surf_107 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_924 = (+surf_107 & -surf_113 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_925 = (+surf_107 & -surf_113 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_926 = (+surf_107 & -surf_113 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_927 = (+surf_107 & -surf_113 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_928 = (+surf_113 & -surf_114 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_929 = (+surf_113 & -surf_114 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_930 = (+surf_113 & -surf_114 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_931 = (+surf_113 & -surf_114 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_932 = (+surf_114 & -surf_115 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_933 = (+surf_114 & -surf_115 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_934 = (+surf_114 & -surf_115 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_935 = (+surf_114 & -surf_115 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_936 = (+surf_115 & -surf_116 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_937 = (+surf_115 & -surf_116 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_938 = (+surf_115 & -surf_116 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_939 = (+surf_115 & -surf_116 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_940 = (+surf_116 & -surf_117 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_941 = (+surf_116 & -surf_117 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_942 = (+surf_116 & -surf_117 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_943 = (+surf_116 & -surf_117 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_944 = (+surf_117 & -surf_118 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_945 = (+surf_117 & -surf_118 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_946 = (+surf_117 & -surf_118 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_947 = (+surf_117 & -surf_118 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_948 = (+surf_118 & -surf_119 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_949 = (+surf_118 & -surf_119 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_950 = (+surf_118 & -surf_119 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_951 = (+surf_118 & -surf_119 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_952 = (+surf_119 & -surf_120 & +surf_85 & -
                  surf_95 & -surf_83 & +surf_93)
    region_953 = (+surf_119 & -surf_120 & -surf_84 & +
                  surf_94 & -surf_83 & +surf_93)
    region_954 = (+surf_119 & -surf_120 & +surf_85 & -
                  surf_95 & +surf_82 & -surf_92)
    region_955 = (+surf_119 & -surf_120 & -surf_84 & +
                  surf_94 & +surf_82 & -surf_92)
    region_990 = (+surf_303 & -surf_101 & +surf_92 & -
                  surf_93 & +surf_84 & -surf_86 & +surf_1)
    region_991 = (+surf_303 & -surf_101 & +surf_92 & -
                  surf_93 & +surf_87 & -surf_85 & +surf_1)
    region_992 = (+surf_303 & -surf_101 & +surf_92 & -
                  surf_93 & +surf_86 & -surf_87 & +surf_1)
    region_993 = (+surf_303 & -surf_101 & +surf_92 & -
                  surf_93 & +surf_85 & -surf_95 & +surf_1)
    region_994 = (+surf_303 & -surf_101 & +surf_92 & -
                  surf_93 & +surf_94 & -surf_84 & +surf_1)

    # cells
    cell_1 = openmc.Cell(cell_id=1, region=region_1, fill=mat_3)
    cell_2 = openmc.Cell(cell_id=2, region=region_2, fill=mat_2)
    cell_3 = openmc.Cell(cell_id=3, region=region_3, fill=mat_1)
    cell_4 = openmc.Cell(cell_id=4, region=region_4, fill=mat_1)
    cell_5 = openmc.Cell(cell_id=5, region=region_5, fill=mat_2)
    cell_6 = openmc.Cell(cell_id=6, region=region_6, fill=mat_1)
    cell_7 = openmc.Cell(cell_id=7, region=region_7, fill=mat_1)
    cell_8 = openmc.Cell(cell_id=8, region=region_8, fill=mat_1)
    cell_9 = openmc.Cell(cell_id=9, region=region_9, fill=mat_1)
    cell_10 = openmc.Cell(cell_id=10, region=region_10, fill=None)
    cell_11 = openmc.Cell(cell_id=11, region=region_11, fill=mat_1)
    cell_12 = openmc.Cell(cell_id=12, region=region_12, fill=mat_2)
    cell_13 = openmc.Cell(cell_id=13, region=region_13, fill=mat_1)
    cell_14 = openmc.Cell(cell_id=14, region=region_14, fill=mat_2)
    cell_15 = openmc.Cell(cell_id=15, region=region_15, fill=mat_1)
    cell_16 = openmc.Cell(cell_id=16, region=region_16, fill=mat_4)
    cell_17 = openmc.Cell(cell_id=17, region=region_17, fill=mat_4)
    cell_18 = openmc.Cell(cell_id=18, region=region_18, fill=mat_4)
    cell_19 = openmc.Cell(cell_id=19, region=region_19, fill=mat_4)
    cell_20 = openmc.Cell(cell_id=20, region=region_20, fill=mat_4)
    cell_21 = openmc.Cell(cell_id=21, region=region_21, fill=mat_4)
    cell_22 = openmc.Cell(cell_id=22, region=region_22, fill=mat_1)
    cell_23 = openmc.Cell(cell_id=23, region=region_23, fill=mat_2)
    cell_24 = openmc.Cell(cell_id=24, region=region_24, fill=mat_2)
    cell_25 = openmc.Cell(cell_id=25, region=region_25, fill=mat_2)
    cell_26 = openmc.Cell(cell_id=26, region=region_26, fill=mat_1)
    cell_27 = openmc.Cell(cell_id=27, region=region_27, fill=mat_2)
    cell_28 = openmc.Cell(cell_id=28, region=region_28, fill=mat_2)
    cell_29 = openmc.Cell(cell_id=29, region=region_29, fill=mat_1)
    cell_30 = openmc.Cell(cell_id=30, region=region_30, fill=mat_2)
    cell_31 = openmc.Cell(cell_id=31, region=region_31, fill=mat_2)
    cell_32 = openmc.Cell(cell_id=32, region=region_32, fill=mat_1)
    cell_33 = openmc.Cell(cell_id=33, region=region_33, fill=mat_1)
    cell_34 = openmc.Cell(cell_id=34, region=region_34, fill=mat_1)
    cell_35 = openmc.Cell(cell_id=35, region=region_35, fill=mat_2)
    cell_36 = openmc.Cell(cell_id=36, region=region_36, fill=mat_2)
    cell_37 = openmc.Cell(cell_id=37, region=region_37, fill=mat_2)
    cell_38 = openmc.Cell(cell_id=38, region=region_38, fill=mat_1)
    cell_39 = openmc.Cell(cell_id=39, region=region_39, fill=mat_2)
    cell_40 = openmc.Cell(cell_id=40, region=region_40, fill=mat_2)
    cell_41 = openmc.Cell(cell_id=41, region=region_41, fill=mat_1)
    cell_42 = openmc.Cell(cell_id=42, region=region_42, fill=mat_2)
    cell_43 = openmc.Cell(cell_id=43, region=region_43, fill=mat_2)
    cell_44 = openmc.Cell(cell_id=44, region=region_44, fill=mat_1)
    cell_45 = openmc.Cell(cell_id=45, region=region_45, fill=mat_1)
    cell_46 = openmc.Cell(cell_id=46, region=region_46, fill=mat_4)
    cell_47 = openmc.Cell(cell_id=47, region=region_47, fill=mat_4)
    cell_48 = openmc.Cell(cell_id=48, region=region_48, fill=mat_4)
    cell_49 = openmc.Cell(cell_id=49, region=region_49, fill=mat_4)
    cell_50 = openmc.Cell(cell_id=50, region=region_50, fill=mat_1)
    cell_51 = openmc.Cell(cell_id=51, region=region_51, fill=mat_4)
    cell_52 = openmc.Cell(cell_id=52, region=region_52, fill=mat_4)
    cell_53 = openmc.Cell(cell_id=53, region=region_53, fill=mat_4)
    cell_54 = openmc.Cell(cell_id=54, region=region_54, fill=mat_4)
    cell_55 = openmc.Cell(cell_id=55, region=region_55, fill=mat_4)
    cell_56 = openmc.Cell(cell_id=56, region=region_56, fill=mat_4)
    cell_57 = openmc.Cell(cell_id=57, region=region_57, fill=mat_4)
    cell_58 = openmc.Cell(cell_id=58, region=region_58, fill=mat_4)
    cell_59 = openmc.Cell(cell_id=59, region=region_59, fill=mat_4)
    cell_60 = openmc.Cell(cell_id=60, region=region_60, fill=mat_4)
    cell_61 = openmc.Cell(cell_id=61, region=region_61, fill=mat_4)
    cell_62 = openmc.Cell(cell_id=62, region=region_62, fill=mat_4)
    cell_63 = openmc.Cell(cell_id=63, region=region_63, fill=mat_4)
    cell_64 = openmc.Cell(cell_id=64, region=region_64, fill=mat_4)
    cell_101 = openmc.Cell(cell_id=101, region=region_101, fill=mat_5)
    cell_102 = openmc.Cell(cell_id=102, region=region_102, fill=mat_5)
    cell_103 = openmc.Cell(cell_id=103, region=region_103, fill=mat_5)
    cell_104 = openmc.Cell(cell_id=104, region=region_104, fill=mat_5)
    cell_105 = openmc.Cell(cell_id=105, region=region_105, fill=mat_5)
    cell_106 = openmc.Cell(cell_id=106, region=region_106, fill=mat_5)
    cell_107 = openmc.Cell(cell_id=107, region=region_107, fill=mat_5)
    cell_108 = openmc.Cell(cell_id=108, region=region_108, fill=mat_5)
    cell_109 = openmc.Cell(cell_id=109, region=region_109, fill=mat_5)
    cell_110 = openmc.Cell(cell_id=110, region=region_110, fill=mat_5)
    cell_111 = openmc.Cell(cell_id=111, region=region_111, fill=mat_5)
    cell_112 = openmc.Cell(cell_id=112, region=region_112, fill=mat_5)
    cell_113 = openmc.Cell(cell_id=113, region=region_113, fill=mat_5)
    cell_114 = openmc.Cell(cell_id=114, region=region_114, fill=mat_5)
    cell_115 = openmc.Cell(cell_id=115, region=region_115, fill=mat_5)
    cell_116 = openmc.Cell(cell_id=116, region=region_116, fill=mat_5)
    cell_117 = openmc.Cell(cell_id=117, region=region_117, fill=mat_5)
    cell_118 = openmc.Cell(cell_id=118, region=region_118, fill=mat_5)
    cell_119 = openmc.Cell(cell_id=119, region=region_119, fill=mat_5)
    cell_120 = openmc.Cell(cell_id=120, region=region_120, fill=mat_5)
    cell_121 = openmc.Cell(cell_id=121, region=region_121, fill=mat_5)
    cell_122 = openmc.Cell(cell_id=122, region=region_122, fill=mat_5)
    cell_124 = openmc.Cell(cell_id=124, region=region_124, fill=mat_5)
    cell_125 = openmc.Cell(cell_id=125, region=region_125, fill=mat_5)
    cell_126 = openmc.Cell(cell_id=126, region=region_126, fill=mat_5)
    cell_127 = openmc.Cell(cell_id=127, region=region_127, fill=mat_5)
    cell_128 = openmc.Cell(cell_id=128, region=region_128, fill=mat_5)
    cell_129 = openmc.Cell(cell_id=129, region=region_129, fill=mat_5)
    cell_130 = openmc.Cell(cell_id=130, region=region_130, fill=mat_5)
    cell_131 = openmc.Cell(cell_id=131, region=region_131, fill=mat_5)
    cell_132 = openmc.Cell(cell_id=132, region=region_132, fill=mat_5)
    cell_133 = openmc.Cell(cell_id=133, region=region_133, fill=mat_5)
    cell_136 = openmc.Cell(cell_id=136, region=region_136, fill=mat_5)
    cell_137 = openmc.Cell(cell_id=137, region=region_137, fill=mat_5)
    cell_138 = openmc.Cell(cell_id=138, region=region_138, fill=mat_5)
    cell_139 = openmc.Cell(cell_id=139, region=region_139, fill=mat_5)
    cell_140 = openmc.Cell(cell_id=140, region=region_140, fill=mat_5)
    cell_141 = openmc.Cell(cell_id=141, region=region_141, fill=mat_5)
    cell_142 = openmc.Cell(cell_id=142, region=region_142, fill=mat_5)
    cell_143 = openmc.Cell(cell_id=143, region=region_143, fill=mat_5)
    cell_144 = openmc.Cell(cell_id=144, region=region_144, fill=mat_5)
    cell_145 = openmc.Cell(cell_id=145, region=region_145, fill=mat_5)
    cell_147 = openmc.Cell(cell_id=147, region=region_147, fill=mat_5)
    cell_148 = openmc.Cell(cell_id=148, region=region_148, fill=mat_5)
    cell_149 = openmc.Cell(cell_id=149, region=region_149, fill=mat_5)
    cell_150 = openmc.Cell(cell_id=150, region=region_150, fill=mat_5)
    cell_151 = openmc.Cell(cell_id=151, region=region_151, fill=mat_5)
    cell_152 = openmc.Cell(cell_id=152, region=region_152, fill=mat_5)
    cell_153 = openmc.Cell(cell_id=153, region=region_153, fill=mat_5)
    cell_154 = openmc.Cell(cell_id=154, region=region_154, fill=mat_5)
    cell_155 = openmc.Cell(cell_id=155, region=region_155, fill=mat_5)
    cell_156 = openmc.Cell(cell_id=156, region=region_156, fill=mat_5)
    cell_157 = openmc.Cell(cell_id=157, region=region_157, fill=mat_5)
    cell_158 = openmc.Cell(cell_id=158, region=region_158, fill=mat_5)
    cell_159 = openmc.Cell(cell_id=159, region=region_159, fill=mat_5)
    cell_160 = openmc.Cell(cell_id=160, region=region_160, fill=mat_5)
    cell_161 = openmc.Cell(cell_id=161, region=region_161, fill=mat_5)
    cell_162 = openmc.Cell(cell_id=162, region=region_162, fill=mat_5)
    cell_163 = openmc.Cell(cell_id=163, region=region_163, fill=mat_5)
    cell_164 = openmc.Cell(cell_id=164, region=region_164, fill=mat_5)
    cell_165 = openmc.Cell(cell_id=165, region=region_165, fill=mat_5)
    cell_166 = openmc.Cell(cell_id=166, region=region_166, fill=mat_5)
    cell_167 = openmc.Cell(cell_id=167, region=region_167, fill=mat_5)
    cell_168 = openmc.Cell(cell_id=168, region=region_168, fill=mat_5)
    cell_228 = openmc.Cell(cell_id=228, region=region_228, fill=mat_5)
    cell_229 = openmc.Cell(cell_id=229, region=region_229, fill=mat_5)
    cell_230 = openmc.Cell(cell_id=230, region=region_230, fill=mat_5)
    cell_231 = openmc.Cell(cell_id=231, region=region_231, fill=mat_5)
    cell_232 = openmc.Cell(cell_id=232, region=region_232, fill=mat_5)
    cell_233 = openmc.Cell(cell_id=233, region=region_233, fill=mat_5)
    cell_234 = openmc.Cell(cell_id=234, region=region_234, fill=mat_5)
    cell_235 = openmc.Cell(cell_id=235, region=region_235, fill=mat_5)
    cell_236 = openmc.Cell(cell_id=236, region=region_236, fill=mat_5)
    cell_237 = openmc.Cell(cell_id=237, region=region_237, fill=mat_5)
    cell_240 = openmc.Cell(cell_id=240, region=region_240, fill=mat_5)
    cell_241 = openmc.Cell(cell_id=241, region=region_241, fill=mat_5)
    cell_242 = openmc.Cell(cell_id=242, region=region_242, fill=mat_5)
    cell_243 = openmc.Cell(cell_id=243, region=region_243, fill=mat_5)
    cell_244 = openmc.Cell(cell_id=244, region=region_244, fill=mat_5)
    cell_245 = openmc.Cell(cell_id=245, region=region_245, fill=mat_5)
    cell_246 = openmc.Cell(cell_id=246, region=region_246, fill=mat_5)
    cell_247 = openmc.Cell(cell_id=247, region=region_247, fill=mat_5)
    cell_248 = openmc.Cell(cell_id=248, region=region_248, fill=mat_5)
    cell_249 = openmc.Cell(cell_id=249, region=region_249, fill=mat_5)
    cell_250 = openmc.Cell(cell_id=250, region=region_250, fill=mat_5)
    cell_251 = openmc.Cell(cell_id=251, region=region_251, fill=mat_5)
    cell_252 = openmc.Cell(cell_id=252, region=region_252, fill=mat_5)
    cell_253 = openmc.Cell(cell_id=253, region=region_253, fill=mat_5)
    cell_254 = openmc.Cell(cell_id=254, region=region_254, fill=mat_5)
    cell_255 = openmc.Cell(cell_id=255, region=region_255, fill=mat_5)
    cell_256 = openmc.Cell(cell_id=256, region=region_256, fill=mat_5)
    cell_257 = openmc.Cell(cell_id=257, region=region_257, fill=mat_5)
    cell_258 = openmc.Cell(cell_id=258, region=region_258, fill=mat_5)
    cell_259 = openmc.Cell(cell_id=259, region=region_259, fill=mat_5)
    cell_260 = openmc.Cell(cell_id=260, region=region_260, fill=mat_5)
    cell_261 = openmc.Cell(cell_id=261, region=region_261, fill=mat_5)
    cell_263 = openmc.Cell(cell_id=263, region=region_263, fill=mat_5)
    cell_264 = openmc.Cell(cell_id=264, region=region_264, fill=mat_5)
    cell_265 = openmc.Cell(cell_id=265, region=region_265, fill=mat_5)
    cell_266 = openmc.Cell(cell_id=266, region=region_266, fill=mat_5)
    cell_267 = openmc.Cell(cell_id=267, region=region_267, fill=mat_5)
    cell_268 = openmc.Cell(cell_id=268, region=region_268, fill=mat_5)
    cell_269 = openmc.Cell(cell_id=269, region=region_269, fill=mat_5)
    cell_270 = openmc.Cell(cell_id=270, region=region_270, fill=mat_5)
    cell_271 = openmc.Cell(cell_id=271, region=region_271, fill=mat_5)
    cell_272 = openmc.Cell(cell_id=272, region=region_272, fill=mat_5)
    cell_274 = openmc.Cell(cell_id=274, region=region_274, fill=mat_5)
    cell_275 = openmc.Cell(cell_id=275, region=region_275, fill=mat_5)
    cell_276 = openmc.Cell(cell_id=276, region=region_276, fill=mat_5)
    cell_277 = openmc.Cell(cell_id=277, region=region_277, fill=mat_5)
    cell_278 = openmc.Cell(cell_id=278, region=region_278, fill=mat_5)
    cell_279 = openmc.Cell(cell_id=279, region=region_279, fill=mat_5)
    cell_280 = openmc.Cell(cell_id=280, region=region_280, fill=mat_5)
    cell_281 = openmc.Cell(cell_id=281, region=region_281, fill=mat_5)
    cell_282 = openmc.Cell(cell_id=282, region=region_282, fill=mat_5)
    cell_283 = openmc.Cell(cell_id=283, region=region_283, fill=mat_5)
    cell_284 = openmc.Cell(cell_id=284, region=region_284, fill=mat_5)
    cell_285 = openmc.Cell(cell_id=285, region=region_285, fill=mat_5)
    cell_286 = openmc.Cell(cell_id=286, region=region_286, fill=mat_5)
    cell_287 = openmc.Cell(cell_id=287, region=region_287, fill=mat_5)
    cell_288 = openmc.Cell(cell_id=288, region=region_288, fill=mat_5)
    cell_289 = openmc.Cell(cell_id=289, region=region_289, fill=mat_5)
    cell_290 = openmc.Cell(cell_id=290, region=region_290, fill=mat_5)
    cell_291 = openmc.Cell(cell_id=291, region=region_291, fill=mat_5)
    cell_292 = openmc.Cell(cell_id=292, region=region_292, fill=mat_5)
    cell_293 = openmc.Cell(cell_id=293, region=region_293, fill=mat_5)
    cell_294 = openmc.Cell(cell_id=294, region=region_294, fill=mat_5)
    cell_295 = openmc.Cell(cell_id=295, region=region_295, fill=mat_5)
    cell_297 = openmc.Cell(cell_id=297, region=region_297, fill=mat_5)
    cell_298 = openmc.Cell(cell_id=298, region=region_298, fill=mat_5)
    cell_299 = openmc.Cell(cell_id=299, region=region_299, fill=mat_5)
    cell_300 = openmc.Cell(cell_id=300, region=region_300, fill=mat_5)
    cell_301 = openmc.Cell(cell_id=301, region=region_301, fill=mat_5)
    cell_302 = openmc.Cell(cell_id=302, region=region_302, fill=mat_5)
    cell_303 = openmc.Cell(cell_id=303, region=region_303, fill=mat_5)
    cell_304 = openmc.Cell(cell_id=304, region=region_304, fill=mat_5)
    cell_305 = openmc.Cell(cell_id=305, region=region_305, fill=mat_5)
    cell_306 = openmc.Cell(cell_id=306, region=region_306, fill=mat_5)
    cell_309 = openmc.Cell(cell_id=309, region=region_309, fill=mat_5)
    cell_310 = openmc.Cell(cell_id=310, region=region_310, fill=mat_5)
    cell_311 = openmc.Cell(cell_id=311, region=region_311, fill=mat_5)
    cell_312 = openmc.Cell(cell_id=312, region=region_312, fill=mat_5)
    cell_313 = openmc.Cell(cell_id=313, region=region_313, fill=mat_5)
    cell_314 = openmc.Cell(cell_id=314, region=region_314, fill=mat_5)
    cell_315 = openmc.Cell(cell_id=315, region=region_315, fill=mat_5)
    cell_316 = openmc.Cell(cell_id=316, region=region_316, fill=mat_5)
    cell_317 = openmc.Cell(cell_id=317, region=region_317, fill=mat_5)
    cell_318 = openmc.Cell(cell_id=318, region=region_318, fill=mat_5)
    cell_404 = openmc.Cell(cell_id=404, region=region_404, fill=mat_6)
    cell_405 = openmc.Cell(cell_id=405, region=region_405, fill=mat_6)
    cell_406 = openmc.Cell(cell_id=406, region=region_406, fill=mat_6)
    cell_407 = openmc.Cell(cell_id=407, region=region_407, fill=mat_6)
    cell_408 = openmc.Cell(cell_id=408, region=region_408, fill=mat_6)
    cell_409 = openmc.Cell(cell_id=409, region=region_409, fill=mat_6)
    cell_411 = openmc.Cell(cell_id=411, region=region_411, fill=mat_6)
    cell_412 = openmc.Cell(cell_id=412, region=region_412, fill=mat_6)
    cell_413 = openmc.Cell(cell_id=413, region=region_413, fill=mat_6)
    cell_414 = openmc.Cell(cell_id=414, region=region_414, fill=mat_6)
    cell_415 = openmc.Cell(cell_id=415, region=region_415, fill=mat_6)
    cell_416 = openmc.Cell(cell_id=416, region=region_416, fill=mat_6)
    cell_417 = openmc.Cell(cell_id=417, region=region_417, fill=mat_6)
    cell_418 = openmc.Cell(cell_id=418, region=region_418, fill=mat_6)
    cell_419 = openmc.Cell(cell_id=419, region=region_419, fill=mat_6)
    cell_421 = openmc.Cell(cell_id=421, region=region_421, fill=mat_6)
    cell_422 = openmc.Cell(cell_id=422, region=region_422, fill=mat_6)
    cell_423 = openmc.Cell(cell_id=423, region=region_423, fill=mat_6)
    cell_424 = openmc.Cell(cell_id=424, region=region_424, fill=mat_6)
    cell_425 = openmc.Cell(cell_id=425, region=region_425, fill=mat_6)
    cell_426 = openmc.Cell(cell_id=426, region=region_426, fill=mat_6)
    cell_427 = openmc.Cell(cell_id=427, region=region_427, fill=mat_6)
    cell_428 = openmc.Cell(cell_id=428, region=region_428, fill=mat_6)
    cell_429 = openmc.Cell(cell_id=429, region=region_429, fill=mat_6)
    cell_431 = openmc.Cell(cell_id=431, region=region_431, fill=mat_6)
    cell_432 = openmc.Cell(cell_id=432, region=region_432, fill=mat_6)
    cell_433 = openmc.Cell(cell_id=433, region=region_433, fill=mat_6)
    cell_434 = openmc.Cell(cell_id=434, region=region_434, fill=mat_6)
    cell_435 = openmc.Cell(cell_id=435, region=region_435, fill=mat_6)
    cell_436 = openmc.Cell(cell_id=436, region=region_436, fill=mat_6)
    cell_437 = openmc.Cell(cell_id=437, region=region_437, fill=mat_6)
    cell_438 = openmc.Cell(cell_id=438, region=region_438, fill=mat_6)
    cell_439 = openmc.Cell(cell_id=439, region=region_439, fill=mat_6)
    cell_441 = openmc.Cell(cell_id=441, region=region_441, fill=mat_6)
    cell_442 = openmc.Cell(cell_id=442, region=region_442, fill=mat_6)
    cell_443 = openmc.Cell(cell_id=443, region=region_443, fill=mat_6)
    cell_444 = openmc.Cell(cell_id=444, region=region_444, fill=mat_6)
    cell_445 = openmc.Cell(cell_id=445, region=region_445, fill=mat_6)
    cell_446 = openmc.Cell(cell_id=446, region=region_446, fill=mat_6)
    cell_447 = openmc.Cell(cell_id=447, region=region_447, fill=mat_6)
    cell_448 = openmc.Cell(cell_id=448, region=region_448, fill=mat_6)
    cell_449 = openmc.Cell(cell_id=449, region=region_449, fill=mat_6)
    cell_504 = openmc.Cell(cell_id=504, region=region_504, fill=mat_6)
    cell_505 = openmc.Cell(cell_id=505, region=region_505, fill=mat_6)
    cell_506 = openmc.Cell(cell_id=506, region=region_506, fill=mat_6)
    cell_507 = openmc.Cell(cell_id=507, region=region_507, fill=mat_6)
    cell_508 = openmc.Cell(cell_id=508, region=region_508, fill=mat_6)
    cell_509 = openmc.Cell(cell_id=509, region=region_509, fill=mat_6)
    cell_511 = openmc.Cell(cell_id=511, region=region_511, fill=mat_6)
    cell_512 = openmc.Cell(cell_id=512, region=region_512, fill=mat_6)
    cell_513 = openmc.Cell(cell_id=513, region=region_513, fill=mat_6)
    cell_514 = openmc.Cell(cell_id=514, region=region_514, fill=mat_6)
    cell_515 = openmc.Cell(cell_id=515, region=region_515, fill=mat_6)
    cell_516 = openmc.Cell(cell_id=516, region=region_516, fill=mat_6)
    cell_517 = openmc.Cell(cell_id=517, region=region_517, fill=mat_6)
    cell_518 = openmc.Cell(cell_id=518, region=region_518, fill=mat_6)
    cell_519 = openmc.Cell(cell_id=519, region=region_519, fill=mat_6)
    cell_521 = openmc.Cell(cell_id=521, region=region_521, fill=mat_6)
    cell_522 = openmc.Cell(cell_id=522, region=region_522, fill=mat_6)
    cell_523 = openmc.Cell(cell_id=523, region=region_523, fill=mat_6)
    cell_524 = openmc.Cell(cell_id=524, region=region_524, fill=mat_6)
    cell_525 = openmc.Cell(cell_id=525, region=region_525, fill=mat_6)
    cell_526 = openmc.Cell(cell_id=526, region=region_526, fill=mat_6)
    cell_527 = openmc.Cell(cell_id=527, region=region_527, fill=mat_6)
    cell_528 = openmc.Cell(cell_id=528, region=region_528, fill=mat_6)
    cell_529 = openmc.Cell(cell_id=529, region=region_529, fill=mat_6)
    cell_531 = openmc.Cell(cell_id=531, region=region_531, fill=mat_6)
    cell_532 = openmc.Cell(cell_id=532, region=region_532, fill=mat_6)
    cell_533 = openmc.Cell(cell_id=533, region=region_533, fill=mat_6)
    cell_534 = openmc.Cell(cell_id=534, region=region_534, fill=mat_6)
    cell_535 = openmc.Cell(cell_id=535, region=region_535, fill=mat_6)
    cell_536 = openmc.Cell(cell_id=536, region=region_536, fill=mat_6)
    cell_537 = openmc.Cell(cell_id=537, region=region_537, fill=mat_6)
    cell_538 = openmc.Cell(cell_id=538, region=region_538, fill=mat_6)
    cell_539 = openmc.Cell(cell_id=539, region=region_539, fill=mat_6)
    cell_541 = openmc.Cell(cell_id=541, region=region_541, fill=mat_6)
    cell_542 = openmc.Cell(cell_id=542, region=region_542, fill=mat_6)
    cell_543 = openmc.Cell(cell_id=543, region=region_543, fill=mat_6)
    cell_544 = openmc.Cell(cell_id=544, region=region_544, fill=mat_6)
    cell_545 = openmc.Cell(cell_id=545, region=region_545, fill=mat_6)
    cell_546 = openmc.Cell(cell_id=546, region=region_546, fill=mat_6)
    cell_547 = openmc.Cell(cell_id=547, region=region_547, fill=mat_6)
    cell_548 = openmc.Cell(cell_id=548, region=region_548, fill=mat_6)
    cell_549 = openmc.Cell(cell_id=549, region=region_549, fill=mat_6)
    cell_600 = openmc.Cell(cell_id=600, region=region_600, fill=None)
    cell_601 = openmc.Cell(cell_id=601, region=region_601, fill=mat_4)
    cell_612 = openmc.Cell(cell_id=612, region=region_612, fill=mat_4)
    cell_621 = openmc.Cell(cell_id=621, region=region_621, fill=mat_4)
    cell_631 = openmc.Cell(cell_id=631, region=region_631, fill=mat_4)
    cell_642 = openmc.Cell(cell_id=642, region=region_642, fill=mat_9)
    cell_644 = openmc.Cell(cell_id=644, region=region_644, fill=mat_9)
    cell_646 = openmc.Cell(cell_id=646, region=region_646, fill=mat_9)
    cell_648 = openmc.Cell(cell_id=648, region=region_648, fill=mat_9)
    cell_654 = openmc.Cell(cell_id=654, region=region_654, fill=None)
    cell_655 = openmc.Cell(cell_id=655, region=region_655, fill=mat_8)
    cell_656 = openmc.Cell(cell_id=656, region=region_656, fill=mat_8)
    cell_657 = openmc.Cell(cell_id=657, region=region_657, fill=mat_8)
    cell_658 = openmc.Cell(cell_id=658, region=region_658, fill=mat_8)
    cell_659 = openmc.Cell(cell_id=659, region=region_659, fill=mat_8)
    cell_660 = openmc.Cell(cell_id=660, region=region_660, fill=mat_8)
    cell_661 = openmc.Cell(cell_id=661, region=region_661, fill=mat_8)
    cell_701 = openmc.Cell(cell_id=701, region=region_701, fill=mat_7)
    cell_702 = openmc.Cell(cell_id=702, region=region_702, fill=mat_7)
    cell_703 = openmc.Cell(cell_id=703, region=region_703, fill=mat_7)
    cell_704 = openmc.Cell(cell_id=704, region=region_704, fill=mat_7)
    cell_705 = openmc.Cell(cell_id=705, region=region_705, fill=mat_7)
    cell_706 = openmc.Cell(cell_id=706, region=region_706, fill=mat_7)
    cell_707 = openmc.Cell(cell_id=707, region=region_707, fill=mat_7)
    cell_708 = openmc.Cell(cell_id=708, region=region_708, fill=mat_7)
    cell_709 = openmc.Cell(cell_id=709, region=region_709, fill=mat_7)
    cell_710 = openmc.Cell(cell_id=710, region=region_710, fill=mat_7)
    cell_711 = openmc.Cell(cell_id=711, region=region_711, fill=mat_7)
    cell_712 = openmc.Cell(cell_id=712, region=region_712, fill=mat_7)
    cell_713 = openmc.Cell(cell_id=713, region=region_713, fill=mat_7)
    cell_714 = openmc.Cell(cell_id=714, region=region_714, fill=mat_7)
    cell_715 = openmc.Cell(cell_id=715, region=region_715, fill=mat_7)
    cell_716 = openmc.Cell(cell_id=716, region=region_716, fill=mat_7)
    cell_717 = openmc.Cell(cell_id=717, region=region_717, fill=mat_7)
    cell_718 = openmc.Cell(cell_id=718, region=region_718, fill=mat_7)
    cell_719 = openmc.Cell(cell_id=719, region=region_719, fill=mat_7)
    cell_720 = openmc.Cell(cell_id=720, region=region_720, fill=mat_7)
    cell_721 = openmc.Cell(cell_id=721, region=region_721, fill=mat_7)
    cell_722 = openmc.Cell(cell_id=722, region=region_722, fill=mat_7)
    cell_723 = openmc.Cell(cell_id=723, region=region_723, fill=mat_7)
    cell_724 = openmc.Cell(cell_id=724, region=region_724, fill=mat_7)
    cell_725 = openmc.Cell(cell_id=725, region=region_725, fill=mat_7)
    cell_726 = openmc.Cell(cell_id=726, region=region_726, fill=mat_7)
    cell_727 = openmc.Cell(cell_id=727, region=region_727, fill=mat_7)
    cell_728 = openmc.Cell(cell_id=728, region=region_728, fill=mat_7)
    cell_729 = openmc.Cell(cell_id=729, region=region_729, fill=mat_7)
    cell_730 = openmc.Cell(cell_id=730, region=region_730, fill=mat_7)
    cell_731 = openmc.Cell(cell_id=731, region=region_731, fill=mat_7)
    cell_732 = openmc.Cell(cell_id=732, region=region_732, fill=mat_7)
    cell_733 = openmc.Cell(cell_id=733, region=region_733, fill=mat_7)
    cell_734 = openmc.Cell(cell_id=734, region=region_734, fill=mat_7)
    cell_735 = openmc.Cell(cell_id=735, region=region_735, fill=mat_7)
    cell_736 = openmc.Cell(cell_id=736, region=region_736, fill=mat_7)
    cell_737 = openmc.Cell(cell_id=737, region=region_737, fill=mat_7)
    cell_738 = openmc.Cell(cell_id=738, region=region_738, fill=mat_7)
    cell_739 = openmc.Cell(cell_id=739, region=region_739, fill=mat_7)
    cell_740 = openmc.Cell(cell_id=740, region=region_740, fill=mat_7)
    cell_741 = openmc.Cell(cell_id=741, region=region_741, fill=mat_7)
    cell_742 = openmc.Cell(cell_id=742, region=region_742, fill=mat_7)
    cell_743 = openmc.Cell(cell_id=743, region=region_743, fill=mat_7)
    cell_744 = openmc.Cell(cell_id=744, region=region_744, fill=mat_7)
    cell_745 = openmc.Cell(cell_id=745, region=region_745, fill=mat_7)
    cell_746 = openmc.Cell(cell_id=746, region=region_746, fill=mat_7)
    cell_747 = openmc.Cell(cell_id=747, region=region_747, fill=mat_7)
    cell_748 = openmc.Cell(cell_id=748, region=region_748, fill=mat_7)
    cell_749 = openmc.Cell(cell_id=749, region=region_749, fill=mat_7)
    cell_750 = openmc.Cell(cell_id=750, region=region_750, fill=mat_7)
    cell_751 = openmc.Cell(cell_id=751, region=region_751, fill=mat_7)
    cell_752 = openmc.Cell(cell_id=752, region=region_752, fill=mat_7)
    cell_753 = openmc.Cell(cell_id=753, region=region_753, fill=mat_7)
    cell_754 = openmc.Cell(cell_id=754, region=region_754, fill=mat_7)
    cell_755 = openmc.Cell(cell_id=755, region=region_755, fill=mat_7)
    cell_756 = openmc.Cell(cell_id=756, region=region_756, fill=mat_7)
    cell_757 = openmc.Cell(cell_id=757, region=region_757, fill=mat_7)
    cell_758 = openmc.Cell(cell_id=758, region=region_758, fill=mat_7)
    cell_759 = openmc.Cell(cell_id=759, region=region_759, fill=mat_7)
    cell_760 = openmc.Cell(cell_id=760, region=region_760, fill=mat_7)
    cell_761 = openmc.Cell(cell_id=761, region=region_761, fill=mat_7)
    cell_762 = openmc.Cell(cell_id=762, region=region_762, fill=mat_7)
    cell_763 = openmc.Cell(cell_id=763, region=region_763, fill=mat_7)
    cell_764 = openmc.Cell(cell_id=764, region=region_764, fill=mat_7)
    cell_765 = openmc.Cell(cell_id=765, region=region_765, fill=mat_7)
    cell_766 = openmc.Cell(cell_id=766, region=region_766, fill=mat_7)
    cell_767 = openmc.Cell(cell_id=767, region=region_767, fill=mat_7)
    cell_768 = openmc.Cell(cell_id=768, region=region_768, fill=mat_7)
    cell_769 = openmc.Cell(cell_id=769, region=region_769, fill=mat_7)
    cell_770 = openmc.Cell(cell_id=770, region=region_770, fill=mat_7)
    cell_771 = openmc.Cell(cell_id=771, region=region_771, fill=mat_7)
    cell_772 = openmc.Cell(cell_id=772, region=region_772, fill=mat_7)
    cell_773 = openmc.Cell(cell_id=773, region=region_773, fill=mat_7)
    cell_774 = openmc.Cell(cell_id=774, region=region_774, fill=mat_7)
    cell_775 = openmc.Cell(cell_id=775, region=region_775, fill=mat_7)
    cell_776 = openmc.Cell(cell_id=776, region=region_776, fill=mat_7)
    cell_777 = openmc.Cell(cell_id=777, region=region_777, fill=mat_7)
    cell_778 = openmc.Cell(cell_id=778, region=region_778, fill=mat_7)
    cell_779 = openmc.Cell(cell_id=779, region=region_779, fill=mat_7)
    cell_780 = openmc.Cell(cell_id=780, region=region_780, fill=mat_7)
    cell_781 = openmc.Cell(cell_id=781, region=region_781, fill=mat_7)
    cell_782 = openmc.Cell(cell_id=782, region=region_782, fill=mat_7)
    cell_783 = openmc.Cell(cell_id=783, region=region_783, fill=mat_7)
    cell_784 = openmc.Cell(cell_id=784, region=region_784, fill=mat_7)
    cell_785 = openmc.Cell(cell_id=785, region=region_785, fill=mat_7)
    cell_786 = openmc.Cell(cell_id=786, region=region_786, fill=mat_7)
    cell_787 = openmc.Cell(cell_id=787, region=region_787, fill=mat_7)
    cell_788 = openmc.Cell(cell_id=788, region=region_788, fill=mat_7)
    cell_789 = openmc.Cell(cell_id=789, region=region_789, fill=mat_7)
    cell_790 = openmc.Cell(cell_id=790, region=region_790, fill=mat_7)
    cell_791 = openmc.Cell(cell_id=791, region=region_791, fill=mat_7)
    cell_792 = openmc.Cell(cell_id=792, region=region_792, fill=mat_7)
    cell_793 = openmc.Cell(cell_id=793, region=region_793, fill=mat_7)
    cell_794 = openmc.Cell(cell_id=794, region=region_794, fill=mat_7)
    cell_795 = openmc.Cell(cell_id=795, region=region_795, fill=mat_7)
    cell_796 = openmc.Cell(cell_id=796, region=region_796, fill=mat_7)
    cell_797 = openmc.Cell(cell_id=797, region=region_797, fill=mat_7)
    cell_798 = openmc.Cell(cell_id=798, region=region_798, fill=mat_7)
    cell_801 = openmc.Cell(cell_id=801, region=region_801, fill=mat_7)
    cell_802 = openmc.Cell(cell_id=802, region=region_802, fill=mat_7)
    cell_803 = openmc.Cell(cell_id=803, region=region_803, fill=mat_7)
    cell_804 = openmc.Cell(cell_id=804, region=region_804, fill=mat_7)
    cell_805 = openmc.Cell(cell_id=805, region=region_805, fill=mat_7)
    cell_806 = openmc.Cell(cell_id=806, region=region_806, fill=mat_7)
    cell_807 = openmc.Cell(cell_id=807, region=region_807, fill=mat_7)
    cell_808 = openmc.Cell(cell_id=808, region=region_808, fill=mat_7)
    cell_809 = openmc.Cell(cell_id=809, region=region_809, fill=mat_7)
    cell_810 = openmc.Cell(cell_id=810, region=region_810, fill=mat_7)
    cell_811 = openmc.Cell(cell_id=811, region=region_811, fill=mat_7)
    cell_812 = openmc.Cell(cell_id=812, region=region_812, fill=mat_7)
    cell_813 = openmc.Cell(cell_id=813, region=region_813, fill=mat_7)
    cell_814 = openmc.Cell(cell_id=814, region=region_814, fill=mat_7)
    cell_815 = openmc.Cell(cell_id=815, region=region_815, fill=mat_7)
    cell_816 = openmc.Cell(cell_id=816, region=region_816, fill=mat_7)
    cell_817 = openmc.Cell(cell_id=817, region=region_817, fill=mat_7)
    cell_818 = openmc.Cell(cell_id=818, region=region_818, fill=mat_7)
    cell_819 = openmc.Cell(cell_id=819, region=region_819, fill=mat_7)
    cell_820 = openmc.Cell(cell_id=820, region=region_820, fill=mat_7)
    cell_821 = openmc.Cell(cell_id=821, region=region_821, fill=mat_7)
    cell_822 = openmc.Cell(cell_id=822, region=region_822, fill=mat_7)
    cell_823 = openmc.Cell(cell_id=823, region=region_823, fill=mat_7)
    cell_824 = openmc.Cell(cell_id=824, region=region_824, fill=mat_7)
    cell_825 = openmc.Cell(cell_id=825, region=region_825, fill=mat_7)
    cell_826 = openmc.Cell(cell_id=826, region=region_826, fill=mat_7)
    cell_827 = openmc.Cell(cell_id=827, region=region_827, fill=mat_7)
    cell_828 = openmc.Cell(cell_id=828, region=region_828, fill=mat_7)
    cell_829 = openmc.Cell(cell_id=829, region=region_829, fill=mat_7)
    cell_830 = openmc.Cell(cell_id=830, region=region_830, fill=mat_7)
    cell_831 = openmc.Cell(cell_id=831, region=region_831, fill=mat_7)
    cell_832 = openmc.Cell(cell_id=832, region=region_832, fill=mat_7)
    cell_833 = openmc.Cell(cell_id=833, region=region_833, fill=mat_7)
    cell_834 = openmc.Cell(cell_id=834, region=region_834, fill=mat_7)
    cell_835 = openmc.Cell(cell_id=835, region=region_835, fill=mat_7)
    cell_836 = openmc.Cell(cell_id=836, region=region_836, fill=mat_7)
    cell_837 = openmc.Cell(cell_id=837, region=region_837, fill=mat_7)
    cell_838 = openmc.Cell(cell_id=838, region=region_838, fill=mat_7)
    cell_839 = openmc.Cell(cell_id=839, region=region_839, fill=mat_7)
    cell_840 = openmc.Cell(cell_id=840, region=region_840, fill=mat_7)
    cell_841 = openmc.Cell(cell_id=841, region=region_841, fill=mat_7)
    cell_842 = openmc.Cell(cell_id=842, region=region_842, fill=mat_7)
    cell_843 = openmc.Cell(cell_id=843, region=region_843, fill=mat_7)
    cell_844 = openmc.Cell(cell_id=844, region=region_844, fill=mat_7)
    cell_845 = openmc.Cell(cell_id=845, region=region_845, fill=mat_7)
    cell_846 = openmc.Cell(cell_id=846, region=region_846, fill=mat_7)
    cell_847 = openmc.Cell(cell_id=847, region=region_847, fill=mat_7)
    cell_848 = openmc.Cell(cell_id=848, region=region_848, fill=mat_7)
    cell_849 = openmc.Cell(cell_id=849, region=region_849, fill=mat_7)
    cell_850 = openmc.Cell(cell_id=850, region=region_850, fill=mat_7)
    cell_851 = openmc.Cell(cell_id=851, region=region_851, fill=mat_7)
    cell_852 = openmc.Cell(cell_id=852, region=region_852, fill=mat_7)
    cell_853 = openmc.Cell(cell_id=853, region=region_853, fill=mat_7)
    cell_854 = openmc.Cell(cell_id=854, region=region_854, fill=mat_7)
    cell_855 = openmc.Cell(cell_id=855, region=region_855, fill=mat_7)
    cell_856 = openmc.Cell(cell_id=856, region=region_856, fill=mat_7)
    cell_857 = openmc.Cell(cell_id=857, region=region_857, fill=mat_7)
    cell_858 = openmc.Cell(cell_id=858, region=region_858, fill=mat_7)
    cell_859 = openmc.Cell(cell_id=859, region=region_859, fill=mat_7)
    cell_860 = openmc.Cell(cell_id=860, region=region_860, fill=mat_7)
    cell_861 = openmc.Cell(cell_id=861, region=region_861, fill=mat_7)
    cell_862 = openmc.Cell(cell_id=862, region=region_862, fill=mat_7)
    cell_863 = openmc.Cell(cell_id=863, region=region_863, fill=mat_7)
    cell_864 = openmc.Cell(cell_id=864, region=region_864, fill=mat_7)
    cell_865 = openmc.Cell(cell_id=865, region=region_865, fill=mat_7)
    cell_866 = openmc.Cell(cell_id=866, region=region_866, fill=mat_7)
    cell_867 = openmc.Cell(cell_id=867, region=region_867, fill=mat_7)
    cell_868 = openmc.Cell(cell_id=868, region=region_868, fill=mat_7)
    cell_869 = openmc.Cell(cell_id=869, region=region_869, fill=mat_7)
    cell_870 = openmc.Cell(cell_id=870, region=region_870, fill=mat_7)
    cell_871 = openmc.Cell(cell_id=871, region=region_871, fill=mat_7)
    cell_872 = openmc.Cell(cell_id=872, region=region_872, fill=mat_7)
    cell_873 = openmc.Cell(cell_id=873, region=region_873, fill=mat_7)
    cell_874 = openmc.Cell(cell_id=874, region=region_874, fill=mat_7)
    cell_875 = openmc.Cell(cell_id=875, region=region_875, fill=mat_7)
    cell_876 = openmc.Cell(cell_id=876, region=region_876, fill=mat_7)
    cell_877 = openmc.Cell(cell_id=877, region=region_877, fill=mat_7)
    cell_878 = openmc.Cell(cell_id=878, region=region_878, fill=mat_7)
    cell_879 = openmc.Cell(cell_id=879, region=region_879, fill=mat_7)
    cell_880 = openmc.Cell(cell_id=880, region=region_880, fill=mat_7)
    cell_881 = openmc.Cell(cell_id=881, region=region_881, fill=mat_7)
    cell_882 = openmc.Cell(cell_id=882, region=region_882, fill=mat_7)
    cell_883 = openmc.Cell(cell_id=883, region=region_883, fill=mat_7)
    cell_884 = openmc.Cell(cell_id=884, region=region_884, fill=mat_7)
    cell_885 = openmc.Cell(cell_id=885, region=region_885, fill=mat_7)
    cell_886 = openmc.Cell(cell_id=886, region=region_886, fill=mat_7)
    cell_887 = openmc.Cell(cell_id=887, region=region_887, fill=mat_7)
    cell_888 = openmc.Cell(cell_id=888, region=region_888, fill=mat_7)
    cell_889 = openmc.Cell(cell_id=889, region=region_889, fill=mat_7)
    cell_890 = openmc.Cell(cell_id=890, region=region_890, fill=mat_7)
    cell_891 = openmc.Cell(cell_id=891, region=region_891, fill=mat_7)
    cell_892 = openmc.Cell(cell_id=892, region=region_892, fill=mat_7)
    cell_893 = openmc.Cell(cell_id=893, region=region_893, fill=mat_7)
    cell_894 = openmc.Cell(cell_id=894, region=region_894, fill=mat_7)
    cell_895 = openmc.Cell(cell_id=895, region=region_895, fill=mat_7)
    cell_896 = openmc.Cell(cell_id=896, region=region_896, fill=mat_7)
    cell_897 = openmc.Cell(cell_id=897, region=region_897, fill=mat_7)
    cell_898 = openmc.Cell(cell_id=898, region=region_898, fill=mat_7)
    cell_900 = openmc.Cell(cell_id=900, region=region_900, fill=mat_7)
    cell_901 = openmc.Cell(cell_id=901, region=region_901, fill=mat_7)
    cell_902 = openmc.Cell(cell_id=902, region=region_902, fill=mat_7)
    cell_903 = openmc.Cell(cell_id=903, region=region_903, fill=mat_7)
    cell_904 = openmc.Cell(cell_id=904, region=region_904, fill=mat_7)
    cell_905 = openmc.Cell(cell_id=905, region=region_905, fill=mat_7)
    cell_906 = openmc.Cell(cell_id=906, region=region_906, fill=mat_7)
    cell_907 = openmc.Cell(cell_id=907, region=region_907, fill=mat_7)
    cell_908 = openmc.Cell(cell_id=908, region=region_908, fill=mat_7)
    cell_909 = openmc.Cell(cell_id=909, region=region_909, fill=mat_7)
    cell_910 = openmc.Cell(cell_id=910, region=region_910, fill=mat_7)
    cell_911 = openmc.Cell(cell_id=911, region=region_911, fill=mat_7)
    cell_912 = openmc.Cell(cell_id=912, region=region_912, fill=mat_7)
    cell_913 = openmc.Cell(cell_id=913, region=region_913, fill=mat_7)
    cell_914 = openmc.Cell(cell_id=914, region=region_914, fill=mat_7)
    cell_915 = openmc.Cell(cell_id=915, region=region_915, fill=mat_7)
    cell_916 = openmc.Cell(cell_id=916, region=region_916, fill=mat_7)
    cell_917 = openmc.Cell(cell_id=917, region=region_917, fill=mat_7)
    cell_918 = openmc.Cell(cell_id=918, region=region_918, fill=mat_7)
    cell_919 = openmc.Cell(cell_id=919, region=region_919, fill=mat_7)
    cell_920 = openmc.Cell(cell_id=920, region=region_920, fill=mat_7)
    cell_921 = openmc.Cell(cell_id=921, region=region_921, fill=mat_7)
    cell_922 = openmc.Cell(cell_id=922, region=region_922, fill=mat_7)
    cell_923 = openmc.Cell(cell_id=923, region=region_923, fill=mat_7)
    cell_924 = openmc.Cell(cell_id=924, region=region_924, fill=mat_7)
    cell_925 = openmc.Cell(cell_id=925, region=region_925, fill=mat_7)
    cell_926 = openmc.Cell(cell_id=926, region=region_926, fill=mat_7)
    cell_927 = openmc.Cell(cell_id=927, region=region_927, fill=mat_7)
    cell_928 = openmc.Cell(cell_id=928, region=region_928, fill=mat_7)
    cell_929 = openmc.Cell(cell_id=929, region=region_929, fill=mat_7)
    cell_930 = openmc.Cell(cell_id=930, region=region_930, fill=mat_7)
    cell_931 = openmc.Cell(cell_id=931, region=region_931, fill=mat_7)
    cell_932 = openmc.Cell(cell_id=932, region=region_932, fill=mat_7)
    cell_933 = openmc.Cell(cell_id=933, region=region_933, fill=mat_7)
    cell_934 = openmc.Cell(cell_id=934, region=region_934, fill=mat_7)
    cell_935 = openmc.Cell(cell_id=935, region=region_935, fill=mat_7)
    cell_936 = openmc.Cell(cell_id=936, region=region_936, fill=mat_7)
    cell_937 = openmc.Cell(cell_id=937, region=region_937, fill=mat_7)
    cell_938 = openmc.Cell(cell_id=938, region=region_938, fill=mat_7)
    cell_939 = openmc.Cell(cell_id=939, region=region_939, fill=mat_7)
    cell_940 = openmc.Cell(cell_id=940, region=region_940, fill=mat_7)
    cell_941 = openmc.Cell(cell_id=941, region=region_941, fill=mat_7)
    cell_942 = openmc.Cell(cell_id=942, region=region_942, fill=mat_7)
    cell_943 = openmc.Cell(cell_id=943, region=region_943, fill=mat_7)
    cell_944 = openmc.Cell(cell_id=944, region=region_944, fill=mat_7)
    cell_945 = openmc.Cell(cell_id=945, region=region_945, fill=mat_7)
    cell_946 = openmc.Cell(cell_id=946, region=region_946, fill=mat_7)
    cell_947 = openmc.Cell(cell_id=947, region=region_947, fill=mat_7)
    cell_948 = openmc.Cell(cell_id=948, region=region_948, fill=mat_7)
    cell_949 = openmc.Cell(cell_id=949, region=region_949, fill=mat_7)
    cell_950 = openmc.Cell(cell_id=950, region=region_950, fill=mat_7)
    cell_951 = openmc.Cell(cell_id=951, region=region_951, fill=mat_7)
    cell_952 = openmc.Cell(cell_id=952, region=region_952, fill=mat_7)
    cell_953 = openmc.Cell(cell_id=953, region=region_953, fill=mat_7)
    cell_954 = openmc.Cell(cell_id=954, region=region_954, fill=mat_7)
    cell_955 = openmc.Cell(cell_id=955, region=region_955, fill=mat_7)
    cell_990 = openmc.Cell(cell_id=990, region=region_990, fill=mat_7)
    cell_991 = openmc.Cell(cell_id=991, region=region_991, fill=mat_7)
    cell_992 = openmc.Cell(cell_id=992, region=region_992, fill=mat_4)
    cell_993 = openmc.Cell(cell_id=993, region=region_993, fill=mat_4)
    cell_994 = openmc.Cell(cell_id=994, region=region_994, fill=mat_4)

    # assign materials to foils/detectors
    if args.reaction_rates:
        cell_602 = openmc.Cell(cell_id=602, region=region_602, fill=mat_100)
        cell_611 = openmc.Cell(cell_id=611, region=region_611, fill=mat_100)
        cell_613 = openmc.Cell(cell_id=613, region=region_613, fill=mat_100)
        cell_614 = openmc.Cell(cell_id=614, region=region_614, fill=mat_100)
        cell_622 = openmc.Cell(cell_id=622, region=region_622, fill=mat_100)
        cell_623 = openmc.Cell(cell_id=623, region=region_623, fill=mat_100)
        cell_632 = openmc.Cell(cell_id=632, region=region_632, fill=mat_100)
        cell_633 = openmc.Cell(cell_id=633, region=region_633, fill=mat_100)
        cell_641 = openmc.Cell(cell_id=641, region=region_641, fill=mat_100)
        cell_643 = openmc.Cell(cell_id=643, region=region_643, fill=mat_100)
        cell_645 = openmc.Cell(cell_id=645, region=region_645, fill=mat_100)
        cell_647 = openmc.Cell(cell_id=647, region=region_647, fill=mat_100)
    elif args.heating:
        cell_602 = openmc.Cell(cell_id=602, region=region_602, fill=mat_60)
        cell_611 = openmc.Cell(cell_id=611, region=region_611, fill=mat_61)
        cell_613 = openmc.Cell(cell_id=613, region=region_613, fill=mat_60)
        cell_614 = openmc.Cell(cell_id=614, region=region_614, fill=mat_61)
        cell_622 = openmc.Cell(cell_id=622, region=region_622, fill=mat_60)
        cell_623 = openmc.Cell(cell_id=623, region=region_623, fill=mat_61)
        cell_632 = openmc.Cell(cell_id=632, region=region_632, fill=mat_60)
        cell_633 = openmc.Cell(cell_id=633, region=region_633, fill=mat_61)
        cell_641 = openmc.Cell(cell_id=641, region=region_641, fill=mat_60)
        cell_643 = openmc.Cell(cell_id=643, region=region_643, fill=mat_60)
        cell_645 = openmc.Cell(cell_id=645, region=region_645, fill=mat_60)
        cell_647 = openmc.Cell(cell_id=647, region=region_647, fill=mat_60)

    # create root universe
    universe = openmc.Universe(cells=[cell_1, cell_2, cell_3, cell_4, cell_5,
                                      cell_6, cell_7, cell_8, cell_9, cell_10, cell_11,
                                      cell_12, cell_13, cell_14, cell_15, cell_16, cell_17,
                                      cell_18, cell_19, cell_20, cell_21, cell_22, cell_23,
                                      cell_24, cell_25, cell_26, cell_27, cell_28, cell_29,
                                      cell_30, cell_31, cell_32, cell_33, cell_34, cell_35,
                                      cell_36, cell_37, cell_38, cell_39, cell_40, cell_41,
                                      cell_42, cell_43, cell_44, cell_45, cell_46, cell_47,
                                      cell_48, cell_49, cell_50, cell_51, cell_52, cell_53,
                                      cell_54, cell_55, cell_56, cell_57, cell_58, cell_59,
                                      cell_60, cell_61, cell_62, cell_63, cell_64, cell_101,
                                      cell_102, cell_103, cell_104, cell_105, cell_106,
                                      cell_107, cell_108, cell_109, cell_110, cell_111,
                                      cell_112, cell_113, cell_114, cell_115, cell_116,
                                      cell_117, cell_118, cell_119, cell_120, cell_121,
                                      cell_122, cell_124, cell_125, cell_126, cell_127,
                                      cell_128, cell_129, cell_130, cell_131, cell_132,
                                      cell_133, cell_136, cell_137, cell_138, cell_139,
                                      cell_140, cell_141, cell_142, cell_143, cell_144,
                                      cell_145, cell_147, cell_148, cell_149, cell_150,
                                      cell_151, cell_152, cell_153, cell_154, cell_155,
                                      cell_156, cell_157, cell_158, cell_159, cell_160,
                                      cell_161, cell_162, cell_163, cell_164, cell_165,
                                      cell_166, cell_167, cell_168, cell_228, cell_229,
                                      cell_230, cell_231, cell_232, cell_233, cell_234,
                                      cell_235, cell_236, cell_237, cell_240, cell_241,
                                      cell_242, cell_243, cell_244, cell_245, cell_246,
                                      cell_247, cell_248, cell_249, cell_250, cell_251,
                                      cell_252, cell_253, cell_254, cell_255, cell_256,
                                      cell_257, cell_258, cell_259, cell_260, cell_261,
                                      cell_263, cell_264, cell_265, cell_266, cell_267,
                                      cell_268, cell_269, cell_270, cell_271, cell_272,
                                      cell_274, cell_275, cell_276, cell_277, cell_278,
                                      cell_279, cell_280, cell_281, cell_282, cell_283,
                                      cell_284, cell_285, cell_286, cell_287, cell_288,
                                      cell_289, cell_290, cell_291, cell_292, cell_293,
                                      cell_294, cell_295, cell_297, cell_298, cell_299,
                                      cell_300, cell_301, cell_302, cell_303, cell_304,
                                      cell_305, cell_306, cell_309, cell_310, cell_311,
                                      cell_312, cell_313, cell_314, cell_315, cell_316,
                                      cell_317, cell_318, cell_404, cell_405, cell_406,
                                      cell_407, cell_408, cell_409, cell_411, cell_412,
                                      cell_413, cell_414, cell_415, cell_416, cell_417,
                                      cell_418, cell_419, cell_421, cell_422, cell_423,
                                      cell_424, cell_425, cell_426, cell_427, cell_428,
                                      cell_429, cell_431, cell_432, cell_433, cell_434,
                                      cell_435, cell_436, cell_437, cell_438, cell_439,
                                      cell_441, cell_442, cell_443, cell_444, cell_445,
                                      cell_446, cell_447, cell_448, cell_449, cell_504,
                                      cell_505, cell_506, cell_507, cell_508, cell_509,
                                      cell_511, cell_512, cell_513, cell_514, cell_515,
                                      cell_516, cell_517, cell_518, cell_519, cell_521,
                                      cell_522, cell_523, cell_524, cell_525, cell_526,
                                      cell_527, cell_528, cell_529, cell_531, cell_532,
                                      cell_533, cell_534, cell_535, cell_536, cell_537,
                                      cell_538, cell_539, cell_541, cell_542, cell_543,
                                      cell_544, cell_545, cell_546, cell_547, cell_548,
                                      cell_549, cell_600, cell_601, cell_602, cell_611,
                                      cell_612, cell_613, cell_614, cell_621, cell_622,
                                      cell_623, cell_631, cell_632, cell_633, cell_641,
                                      cell_642, cell_643, cell_644, cell_645, cell_646,
                                      cell_647, cell_648, cell_654, cell_655, cell_656,
                                      cell_657, cell_658, cell_659, cell_660, cell_661,
                                      cell_701, cell_702, cell_703, cell_704, cell_705,
                                      cell_706, cell_707, cell_708, cell_709, cell_710,
                                      cell_711, cell_712, cell_713, cell_714, cell_715,
                                      cell_716, cell_717, cell_718, cell_719, cell_720,
                                      cell_721, cell_722, cell_723, cell_724, cell_725,
                                      cell_726, cell_727, cell_728, cell_729, cell_730,
                                      cell_731, cell_732, cell_733, cell_734, cell_735,
                                      cell_736, cell_737, cell_738, cell_739, cell_740,
                                      cell_741, cell_742, cell_743, cell_744, cell_745,
                                      cell_746, cell_747, cell_748, cell_749, cell_750,
                                      cell_751, cell_752, cell_753, cell_754, cell_755,
                                      cell_756, cell_757, cell_758, cell_759, cell_760,
                                      cell_761, cell_762, cell_763, cell_764, cell_765,
                                      cell_766, cell_767, cell_768, cell_769, cell_770,
                                      cell_771, cell_772, cell_773, cell_774, cell_775,
                                      cell_776, cell_777, cell_778, cell_779, cell_780,
                                      cell_781, cell_782, cell_783, cell_784, cell_785,
                                      cell_786, cell_787, cell_788, cell_789, cell_790,
                                      cell_791, cell_792, cell_793, cell_794, cell_795,
                                      cell_796, cell_797, cell_798, cell_801, cell_802,
                                      cell_803, cell_804, cell_805, cell_806, cell_807,
                                      cell_808, cell_809, cell_810, cell_811, cell_812,
                                      cell_813, cell_814, cell_815, cell_816, cell_817,
                                      cell_818, cell_819, cell_820, cell_821, cell_822,
                                      cell_823, cell_824, cell_825, cell_826, cell_827,
                                      cell_828, cell_829, cell_830, cell_831, cell_832,
                                      cell_833, cell_834, cell_835, cell_836, cell_837,
                                      cell_838, cell_839, cell_840, cell_841, cell_842,
                                      cell_843, cell_844, cell_845, cell_846, cell_847,
                                      cell_848, cell_849, cell_850, cell_851, cell_852,
                                      cell_853, cell_854, cell_855, cell_856, cell_857,
                                      cell_858, cell_859, cell_860, cell_861, cell_862,
                                      cell_863, cell_864, cell_865, cell_866, cell_867,
                                      cell_868, cell_869, cell_870, cell_871, cell_872,
                                      cell_873, cell_874, cell_875, cell_876, cell_877,
                                      cell_878, cell_879, cell_880, cell_881, cell_882,
                                      cell_883, cell_884, cell_885, cell_886, cell_887,
                                      cell_888, cell_889, cell_890, cell_891, cell_892,
                                      cell_893, cell_894, cell_895, cell_896, cell_897,
                                      cell_898, cell_900, cell_901, cell_902, cell_903,
                                      cell_904, cell_905, cell_906, cell_907, cell_908,
                                      cell_909, cell_910, cell_911, cell_912, cell_913,
                                      cell_914, cell_915, cell_916, cell_917, cell_918,
                                      cell_919, cell_920, cell_921, cell_922, cell_923,
                                      cell_924, cell_925, cell_926, cell_927, cell_928,
                                      cell_929, cell_930, cell_931, cell_932, cell_933,
                                      cell_934, cell_935, cell_936, cell_937, cell_938,
                                      cell_939, cell_940, cell_941, cell_942, cell_943,
                                      cell_944, cell_945, cell_946, cell_947, cell_948,
                                      cell_949, cell_950, cell_951, cell_952, cell_953,
                                      cell_954, cell_955, cell_990, cell_991, cell_992,
                                      cell_993, cell_994])

    # create geometry instance
    model.geometry = openmc.Geometry(universe)

    ############################################################################
    # Define Settings

    # source definition
    # fng source
    fng_center = (0, 0, 0)
    fng_uvw = (0., 1., 0)

    source = fng_source(center=fng_center,
                        reference_uvw=fng_uvw)

    # weight windows from wwinps
    ww = openmc.wwinp_to_wws("weight_windows.cadis.wwinp")

    # settings
    settings = openmc.Settings(run_mode='fixed source')
    settings.batches = args.batches
    settings.particles = args.particles
    settings.weight_windows = ww
    settings.source = source
    if args.heating:
        settings.survival_biasing = True
        settings.photon_transport = True
        settings.electron_treatment = 'ttb'
    settings.output = {'tallies': False}

    ############################################################################
    # Specify Tallies

    model.tallies = openmc.Tallies()

    # filters
    # particle filters
    neutron_filter = openmc.ParticleFilter(['neutron'])
    particle_filter = openmc.ParticleFilter(
        ['neutron', 'photon', 'electron', 'positron'])

    # cell filters
    front_cell_filter = openmc.CellFilter(
        [cell_602, cell_613, cell_622, cell_632])
    mid_cell_filter = openmc.CellFilter(
        [cell_611, cell_614, cell_623, cell_633])
    back_cell_filter = openmc.CellFilter(
        [cell_641, cell_643, cell_645, cell_647])
    heatdetector_cell_filter = openmc.CellFilter(
        [cell_611, cell_614, cell_623, cell_633])

    nuclides = ['zr90', 'al27', 'mn55', 'nb93', 'ni58_n2n',
                'ni58_np', 'au197', 'fe56', 'in115']

    # dosimetry tallies from IRDFF-II nuclear data library
    zr90_n2n_acef = irdff.path + "dos-irdff2-4025.acef"
    al27_na_acef = irdff.path + "dos-irdff2-1325.acef"
    mn55_ng_acef = irdff.path + "dos-irdff2-2525_modified.acef"
    nb93_n2n_acef = irdff.path + "dos-irdff2-4125.acef"
    ni58_n2n_acef = irdff.path + "dos-irdff2-2825_modified.acef"
    ni58_np_acef = irdff.path + "dos-irdff2-2825_modified.acef"
    au197_ng_acef = irdff.path + "dos-irdff2-7925_modified.acef"
    fe56_np_acef = irdff.path + "dos-irdff2-2631.acef"
    in115_nn_acef = irdff.path + "dos-irdff2-4931.acef"

    irdff_xs = [zr90_n2n_acef, al27_na_acef, mn55_ng_acef,
                nb93_n2n_acef, ni58_n2n_acef, ni58_np_acef, au197_ng_acef,
                fe56_np_acef, in115_nn_acef]
    reactions = [16, 107, 102, 11016, 16, 103, 102, 103, 11004]

    cell_filter_list = [front_cell_filter, mid_cell_filter, back_cell_filter,
                        front_cell_filter, mid_cell_filter, mid_cell_filter, back_cell_filter,
                        front_cell_filter, mid_cell_filter]

    # define tallies according to simulation type
    if args.reaction_rates:
        for n, r, x, c in zip(nuclides, reactions, irdff_xs, cell_filter_list):
            tally = openmc.Tally(name=f"rr_{n}")
            irdff_xs = irdff.cross_section(x)
            multiplier = openmc.EnergyFunctionFilter.from_tabulated1d(
                irdff_xs[r])
            tally.filters = [c, neutron_filter, multiplier]
            tally.scores = ["flux"]
            model.tallies.extend([tally])
    elif args.heating:
        tally = openmc.Tally(name='nuclear_heating')
        tally.filters = [heatdetector_cell_filter, particle_filter]
        tally.scores = ['heating']
        model.tallies.append(tally)

    # define the folder names for storing the statepoints
    if args.reaction_rates:
        cwd = 'reaction_rates'
    else:
        cwd = 'heating'

    model.settings = settings

    return model.run(cwd=cwd, threads=args.threads)


if __name__ == "__main__":
    main()
