#!/usr/bin/env python3
import numpy as np
import argparse

import openmc
from openmc_fusion_benchmarks import from_irdff as irdff
from openmc_fusion_benchmarks.neutron_sources import fng_source


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batches", type=int, default=100)
    parser.add_argument("-p", "--particles", type=int, default=int(1e9))
    parser.add_argument("-s", "--threads", type=int)
    parser.add_argument("-c", "--cwd", type=str)

    args = parser.parse_args()

    return args


def main():
    """Analysis of FNS-Dogleg Duct experiment"""

    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    # define materials

    # cu
    mat_1 = openmc.Material(material_id=1, name="cu")
    mat_1.set_density("g/cm3", 8.93)
    mat_1.add_nuclide("Cu63", 0.6915, "ao")
    mat_1.add_nuclide("Cu65", 0.3085, "ao")
    # cool_water
    mat_2 = openmc.Material(material_id=2, name="cool_water")
    mat_2.set_density("g/cm3", 1.0)
    mat_2.add_nuclide("H1", 0.6668961214200001, "ao")
    mat_2.add_nuclide("H2", 0.00010387858, "ao")
    mat_2.add_nuclide("O16", 0.3322076598, "ao")
    mat_2.add_nuclide("O17", 0.000126207, "ao")
    mat_2.add_nuclide("O18", 0.0006661332, "ao")
    # watercu_mix
    mat_3 = openmc.Material(material_id=3, name="watercu_mix")
    mat_3.set_density("g/cm3", 5.6)
    mat_3.add_nuclide("H1", 0.27995639280000006, "ao")
    mat_3.add_nuclide("H2", 4.3607200000000005e-05, "ao")
    mat_3.add_nuclide("O16", 0.13966688400000002, "ao")
    mat_3.add_nuclide("O17", 5.3060000000000004e-05, "ao")
    mat_3.add_nuclide("O18", 0.000280056, "ao")
    mat_3.add_nuclide("Cu63", 0.40107, "ao")
    mat_3.add_nuclide("Cu65", 0.17892999999999998, "ao")
    # ss304
    mat_4 = openmc.Material(material_id=4, name="ss304")
    mat_4.set_density("g/cm3", 7.8000000000000025)
    mat_4.add_nuclide("Cr50", 0.007930004583164795, "wo")
    mat_4.add_nuclide("Cr52", 0.1590287884675572, "wo")
    mat_4.add_nuclide("Cr53", 0.01837981495548285, "wo")
    mat_4.add_nuclide("Cr54", 0.004661391993795167, "wo")
    mat_4.add_nuclide("Mn55", 0.01, "wo")
    mat_4.add_nuclide("Fe54", 0.040083463381218674, "wo")
    mat_4.add_nuclide("Fe56", 0.6525008543348935, "wo")
    mat_4.add_nuclide("Fe57", 0.015338617157089918, "wo")
    mat_4.add_nuclide("Fe58", 0.0020770651267980177, "wo")
    mat_4.add_nuclide("Ni58", 0.060477934766930905, "wo")
    mat_4.add_nuclide("Ni60", 0.024098366609306605, "wo")
    mat_4.add_nuclide("Ni61", 0.0010650231767670392, "wo")
    mat_4.add_nuclide("Ni62", 0.0034513400731384525, "wo")
    mat_4.add_nuclide("Ni64", 0.0009073353738570029, "wo")
    # aluminum
    mat_5 = openmc.Material(material_id=5, name="aluminum")
    mat_5.set_density("g/cm3", 2.6899999999999995)
    mat_5.add_nuclide("Al27", 1.0, "ao")
    # air
    mat_6 = openmc.Material(material_id=6, name="air")
    mat_6.set_density("g/cm3", 0.0012249999999999995)
    mat_6.add_nuclide("N14", 0.7771428600000001, "ao")
    mat_6.add_nuclide("N15", 0.00285714, "ao")
    mat_6.add_nuclide("O16", 0.219476532, "ao")
    mat_6.add_nuclide("O17", 8.338e-05, "ao")
    mat_6.add_nuclide("O18", 0.000440088, "ao")
    mat_6.add_nuclide("Ar36", 3.336e-05, "ao")
    mat_6.add_nuclide("Ar38", 6.29e-06, "ao")
    mat_6.add_nuclide("Ar40", 0.00996035, "ao")
    # mortar
    mat_7 = openmc.Material(material_id=7, name="mortar")
    mat_7.set_density("g/cm3", 2.1023627049432485)
    mat_7.add_nuclide("H1", 0.006596172552072, "ao")
    mat_7.add_nuclide("H2", 1.027447928e-06, "ao")
    mat_7.add_nuclide("C12", 0.00073473937834, "ao")
    mat_7.add_nuclide("C13", 8.23062166e-06, "ao")
    mat_7.add_nuclide("O16", 0.037590344208, "ao")
    mat_7.add_nuclide("O17", 1.428072e-05, "ao")
    mat_7.add_nuclide("O18", 7.537507199999999e-05, "ao")
    mat_7.add_nuclide("Na23", 0.00050802, "ao")
    mat_7.add_nuclide("Mg24", 0.0003700749174, "ao")
    mat_7.add_nuclide("Mg25", 4.6967748e-05, "ao")
    mat_7.add_nuclide("Mg26", 5.16973346e-05, "ao")
    mat_7.add_nuclide("Al27", 0.0019408, "ao")
    mat_7.add_nuclide("Si28", 0.0104653017896, "ao")
    mat_7.add_nuclide("Si29", 0.0005313981652, "ao")
    mat_7.add_nuclide("Si30", 0.0003503000452, "ao")
    mat_7.add_nuclide("K39", 0.00031173385087, "ao")
    mat_7.add_nuclide("K40", 3.910959e-08, "ao")
    mat_7.add_nuclide("K41", 2.2497039540000002e-05, "ao")
    mat_7.add_nuclide("Ca40", 0.0036006795629999998, "ao")
    mat_7.add_nuclide("Ca42", 2.4031521e-05, "ao")
    mat_7.add_nuclide("Ca43", 5.014305e-06, "ao")
    mat_7.add_nuclide("Ca44", 7.7480298e-05, "ao")
    mat_7.add_nuclide("Ca46", 1.4857200000000002e-07, "ao")
    mat_7.add_nuclide("Ca48", 6.945740999999999e-06, "ao")
    mat_7.add_nuclide("Fe54", 3.980445e-05, "ao")
    mat_7.add_nuclide("Fe56", 0.00062484474, "ao")
    mat_7.add_nuclide("Fe57", 1.443039e-05, "ao")
    mat_7.add_nuclide("Fe58", 1.92042e-06, "ao")
    # concrete
    mat_8 = openmc.Material(material_id=8, name="concrete")
    mat_8.set_density("g/cm3", 2.3563125502648763)
    mat_8.add_nuclide("H1", 0.005580930690468, "ao")
    mat_8.add_nuclide("H2", 8.693095319999999e-07, "ao")
    mat_8.add_nuclide("C12", 0.0005349079097999999, "ao")
    mat_8.add_nuclide("C13", 5.992090199999999e-06, "ao")
    mat_8.add_nuclide("O16", 0.04304732889, "ao")
    mat_8.add_nuclide("O17", 1.635385e-05, "ao")
    mat_8.add_nuclide("O18", 8.631726e-05, "ao")
    mat_8.add_nuclide("Na23", 0.0007859, "ao")
    mat_8.add_nuclide("Mg24", 0.00030190862400000004, "ao")
    mat_8.add_nuclide("Mg25", 3.8316480000000004e-05, "ao")
    mat_8.add_nuclide("Mg26", 4.2174896000000004e-05, "ao")
    mat_8.add_nuclide("Al27", 0.002637, "ao")
    mat_8.add_nuclide("Si28", 0.013659215608, "ao")
    mat_8.add_nuclide("Si29", 0.000693575996, "ao")
    mat_8.add_nuclide("Si30", 0.000457208396, "ao")
    mat_8.add_nuclide("K39", 0.0004931488328, "ao")
    mat_8.add_nuclide("K40", 6.186959999999999e-08, "ao")
    mat_8.add_nuclide("K41", 3.5589297599999996e-05, "ao")
    mat_8.add_nuclide("Ca40", 0.00248556724, "ao")
    mat_8.add_nuclide("Ca42", 1.658908e-05, "ao")
    mat_8.add_nuclide("Ca43", 3.4614e-06, "ao")
    mat_8.add_nuclide("Ca44", 5.3485039999999996e-05, "ao")
    mat_8.add_nuclide("Ca46", 1.0256e-07, "ao")
    mat_8.add_nuclide("Ca48", 4.7946799999999995e-06, "ao")
    mat_8.add_nuclide("Fe54", 3.4245855e-05, "ao")
    mat_8.add_nuclide("Fe56", 0.000537586686, "ao")
    mat_8.add_nuclide("Fe57", 1.2415221e-05, "ao")
    mat_8.add_nuclide("Fe58", 1.652238e-06, "ao")
    # nb93(0.0)-in115(0.0)-au197(1e-05)-xylene(0.0)-air(0.99999)
    mat_9 = openmc.Material(
        material_id=9, name="nb93(0.0)-in115(0.0)-au197(1e-05)-xylene(0.0)-air(0.99999)"
    )
    mat_9.set_density("g/cm3", 0.0014179877499999988)
    mat_9.add_nuclide("Nb93", 0.0, "ao")
    mat_9.add_nuclide("In113", 0.0, "ao")
    mat_9.add_nuclide("In115", 0.0, "ao")
    mat_9.add_nuclide("Au197", 0.011619924630199533, "ao")
    mat_9.add_nuclide("H1", 0.0, "ao")
    mat_9.add_nuclide("H2", 0.0, "ao")
    mat_9.add_nuclide("C12", 0.0, "ao")
    mat_9.add_nuclide("C13", 0.0, "ao")
    mat_9.add_nuclide("N14", 0.7605074440989131, "ao")
    mat_9.add_nuclide("N15", 0.0027959804441010615, "ao")
    mat_9.add_nuclide("O16", 0.21477844677233904, "ao")
    mat_9.add_nuclide("O17", 8.159517889538017e-05, "ao")
    mat_9.add_nuclide("O18", 0.0004306675352567769, "ao")
    mat_9.add_nuclide("Ar36", 3.26459003112243e-05, "ao")
    mat_9.add_nuclide("Ar38", 6.155357103045588e-06, "ao")
    mat_9.add_nuclide("Ar40", 0.009747140082880782, "ao")
    # iron_ass
    mat_11 = openmc.Material(material_id=11, name="iron_ass")
    mat_11.set_density("g/cm3", 7.831199999999998)
    mat_11.add_nuclide("C12", 0.0006018579292000001, "ao")
    mat_11.add_nuclide("C13", 6.7420708e-06, "ao")
    mat_11.add_nuclide("Si28", 0.00036065, "ao")
    mat_11.add_nuclide("Si29", 1.843e-05, "ao")
    mat_11.add_nuclide("Si30", 1.221e-05, "ao")
    mat_11.add_nuclide("P31", 3.0452e-05, "ao")
    mat_11.add_nuclide("S32", 5.5899e-06, "ao")
    mat_11.add_nuclide("S33", 4.4122e-08, "ao")
    mat_11.add_nuclide("S34", 2.4767e-07, "ao")
    mat_11.add_nuclide("S36", 1.1766e-09, "ao")
    mat_11.add_nuclide("Mn55", 0.00090565, "ao")
    mat_11.add_nuclide("Fe54", 0.0048262, "ao")
    mat_11.add_nuclide("Fe56", 0.07632, "ao")
    mat_11.add_nuclide("Fe57", 0.0018306, "ao")
    mat_11.add_nuclide("Fe58", 0.00023299, "ao")
    # iron_aux
    mat_12 = openmc.Material(material_id=12, name="iron_aux")
    mat_12.set_density("g/cm3", 6.459999999999999)
    mat_12.add_nuclide("H1", 0.01631145925764, "ao")
    mat_12.add_nuclide("H2", 2.5407423599999996e-06, "ao")
    mat_12.add_nuclide("C12", 0.008555065329799999, "ao")
    mat_12.add_nuclide("C13", 9.58346702e-05, "ao")
    mat_12.add_nuclide("Si28", 0.00021947, "ao")
    mat_12.add_nuclide("Si29", 1.1113e-05, "ao")
    mat_12.add_nuclide("Si30", 7.3768e-06, "ao")
    mat_12.add_nuclide("P31", 1.9505e-05, "ao")
    mat_12.add_nuclide("S32", 6.7134e-06, "ao")
    mat_12.add_nuclide("S33", 5.2989e-08, "ao")
    mat_12.add_nuclide("S34", 2.9744e-07, "ao")
    mat_12.add_nuclide("S36", 1.413e-09, "ao")
    mat_12.add_nuclide("Mn55", 0.00056358, "ao")
    mat_12.add_nuclide("Fe54", 0.0038754, "ao")
    mat_12.add_nuclide("Fe56", 0.061285, "ao")
    mat_12.add_nuclide("Fe57", 0.00147, "ao")
    mat_12.add_nuclide("Fe58", 0.00018709, "ao")
    # lat_water
    mat_21 = openmc.Material(material_id=21, name="lat_water")
    mat_21.set_density("g/cm3", 0.1)
    mat_21.add_nuclide("H1", 0.6668961214200001, "ao")
    mat_21.add_nuclide("H2", 0.00010387858, "ao")
    mat_21.add_nuclide("O16", 0.3322076598, "ao")
    mat_21.add_nuclide("O17", 0.000126207, "ao")
    mat_21.add_nuclide("O18", 0.0006661332, "ao")
    # bot_water
    mat_22 = openmc.Material(material_id=22, name="bot_water")
    mat_22.set_density("g/cm3", 0.5)
    mat_22.add_nuclide("H1", 0.6668961214200001, "ao")
    mat_22.add_nuclide("H2", 0.00010387858, "ao")
    mat_22.add_nuclide("O16", 0.3322076598, "ao")
    mat_22.add_nuclide("O17", 0.000126207, "ao")
    mat_22.add_nuclide("O18", 0.0006661332, "ao")
    # b10
    mat_31 = openmc.Material(material_id=31, name="b10")
    mat_31.set_density("g/cm3", 0.9999999999999999)
    mat_31.add_nuclide("B10", 1.0, "ao")
    # al27
    mat_32 = openmc.Material(material_id=32, name="al27")
    mat_32.set_density("g/cm3", 2.6999999999999997)
    mat_32.add_nuclide("Al27", 1.0, "ao")
    # nb93
    mat_33 = openmc.Material(material_id=33, name="nb93")
    mat_33.set_density("g/cm3", 8.57)
    mat_33.add_nuclide("Nb93", 1.0, "ao")
    # in115
    mat_34 = openmc.Material(material_id=34, name="in115")
    mat_34.set_density("g/cm3", 7.309999999999999)
    mat_34.add_nuclide("In113", 0.04281, "ao")
    mat_34.add_nuclide("In115", 0.95719, "ao")
    # au197
    mat_35 = openmc.Material(material_id=35, name="au197")
    mat_35.set_density("g/cm3", 19.299999999999997)
    mat_35.add_nuclide("Au197", 1.0, "ao")
    # u235
    mat_36 = openmc.Material(material_id=36, name="u235")
    mat_36.set_density("g/cm3", 0.9999999999999999)
    mat_36.add_nuclide("U235", 1.0, "ao")
    # u238
    mat_37 = openmc.Material(material_id=37, name="u238")
    mat_37.set_density("g/cm3", 0.9999999999999999)
    mat_37.add_nuclide("U238", 1.0, "ao")
    # cd_cover
    mat_38 = openmc.Material(material_id=38, name="cd_cover")
    mat_38.set_density("g/cm3", 0.9999999999999999)
    mat_38.add_nuclide("Cd106", 0.0125, "ao")
    mat_38.add_nuclide("Cd108", 0.0089, "ao")
    mat_38.add_nuclide("Cd110", 0.1249, "ao")
    mat_38.add_nuclide("Cd111", 0.128, "ao")
    mat_38.add_nuclide("Cd112", 0.2413, "ao")
    mat_38.add_nuclide("Cd113", 0.1222, "ao")
    mat_38.add_nuclide("Cd114", 0.2873, "ao")
    mat_38.add_nuclide("Cd116", 0.0749, "ao")
    # xylene
    mat_40 = openmc.Material(material_id=40, name="xylene")
    mat_40.set_density("g/cm3", 0.866)
    mat_40.add_nuclide("H1", 9.9984426, "ao")
    mat_40.add_nuclide("H2", 0.0015574, "ao")
    mat_40.add_nuclide("C12", 7.911376, "ao")
    mat_40.add_nuclide("C13", 0.088624, "ao")
    # ss304_acc1
    mat_41 = openmc.Material(material_id=41, name="ss304_acc1")
    mat_41.set_density("g/cm3", 2.600000000000001)
    mat_41.add_nuclide("Cr50", 0.007930004583164795, "wo")
    mat_41.add_nuclide("Cr52", 0.1590287884675572, "wo")
    mat_41.add_nuclide("Cr53", 0.01837981495548285, "wo")
    mat_41.add_nuclide("Cr54", 0.004661391993795167, "wo")
    mat_41.add_nuclide("Mn55", 0.01, "wo")
    mat_41.add_nuclide("Fe54", 0.040083463381218674, "wo")
    mat_41.add_nuclide("Fe56", 0.6525008543348935, "wo")
    mat_41.add_nuclide("Fe57", 0.015338617157089918, "wo")
    mat_41.add_nuclide("Fe58", 0.0020770651267980177, "wo")
    mat_41.add_nuclide("Ni58", 0.060477934766930905, "wo")
    mat_41.add_nuclide("Ni60", 0.024098366609306605, "wo")
    mat_41.add_nuclide("Ni61", 0.0010650231767670392, "wo")
    mat_41.add_nuclide("Ni62", 0.0034513400731384525, "wo")
    mat_41.add_nuclide("Ni64", 0.0009073353738570029, "wo")
    # ss304_acc2
    mat_42 = openmc.Material(material_id=42, name="ss304_acc2")
    mat_42.set_density("g/cm3", 1.56)
    mat_42.add_nuclide("Cr50", 0.007930004583164795, "wo")
    mat_42.add_nuclide("Cr52", 0.1590287884675572, "wo")
    mat_42.add_nuclide("Cr53", 0.01837981495548285, "wo")
    mat_42.add_nuclide("Cr54", 0.004661391993795167, "wo")
    mat_42.add_nuclide("Mn55", 0.01, "wo")
    mat_42.add_nuclide("Fe54", 0.040083463381218674, "wo")
    mat_42.add_nuclide("Fe56", 0.6525008543348935, "wo")
    mat_42.add_nuclide("Fe57", 0.015338617157089918, "wo")
    mat_42.add_nuclide("Fe58", 0.0020770651267980177, "wo")
    mat_42.add_nuclide("Ni58", 0.060477934766930905, "wo")
    mat_42.add_nuclide("Ni60", 0.024098366609306605, "wo")
    mat_42.add_nuclide("Ni61", 0.0010650231767670392, "wo")
    mat_42.add_nuclide("Ni62", 0.0034513400731384525, "wo")
    mat_42.add_nuclide("Ni64", 0.0009073353738570029, "wo")

    ############################################################################
    # Build Geometry

    #  surfaces
    surf_3 = openmc.Sphere(
        surface_id=3, x0=0.0, y0=-22.86, z0=0.0, r=22.8599, name="so_target1"
    )
    surf_6 = openmc.Sphere(
        surface_id=6, x0=0.0, y0=-22.86, z0=0.0, r=22.94, name="so_target2"
    )
    surf_18 = openmc.Quadric(
        surface_id=18,
        a=0.8768777442929967,
        b=0.01085234050855799,
        c=0.9438650424007774,
        d=0.49999999999999994,
        e=-0.0,
        f=0.0,
        g=12.22338849866847,
        h=-2.4647971791202616,
        j=-0.0,
        k=-103.86555859873886,
        name="ky_target6",
    )
    surf_21 = openmc.Plane(
        surface_id=21,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-12.081064388968102,
        name="py_target7",
    )
    surf_9 = openmc.Sphere(
        surface_id=9, x0=0.0, y0=-22.86, z0=0.0, r=23.04, name="so_target3"
    )
    surf_76 = openmc.YCone(
        surface_id=76, x0=0.0, y0=21.88, z0=0.0, r2=0.00105872, name="ky_collimator51"
    )
    surf_77 = openmc.YCone(
        surface_id=77, x0=0.0, y0=34.17, z0=0.0, r2=0.00105872, name="ky_collimator52"
    )
    surf_78 = openmc.YPlane(surface_id=78, y0=-1.17, name="py_collimator53")
    surf_79 = openmc.YPlane(surface_id=79, y0=-15.0, name="py_collimator54")
    surf_15 = openmc.Quadric(
        surface_id=15,
        a=0.8768777442929967,
        b=0.01085234050855799,
        c=0.9438650424007774,
        d=0.49999999999999994,
        e=-0.0,
        f=0.0,
        g=12.211126194740915,
        h=-2.4190336378447572,
        j=-0.0,
        k=-101.53579448631787,
        name="ky_target5",
    )
    surf_57 = openmc.Plane(
        surface_id=57,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-18.221064388968102,
        name="py_supp31",
    )
    surf_63 = openmc.Quadric(
        surface_id=63,
        a=0.9330127018922217,
        b=0.06698729810778066,
        c=1.0000000000000024,
        d=0.5000000000000006,
        e=-0.0,
        f=0.0,
        g=11.430000000000012,
        h=3.0626592694877317,
        j=0.0,
        k=-45.993804549755424,
        name="cy_supp33",
    )
    surf_60 = openmc.Plane(
        surface_id=60,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-24.221064388968102,
        name="py_supp32",
    )
    surf_12 = openmc.Sphere(
        surface_id=12, x0=0.0, y0=-22.86, z0=0.0, r=23.12, name="so_target4"
    )
    surf_27 = openmc.Quadric(
        surface_id=27,
        a=0.9330127018922217,
        b=0.06698729810778066,
        c=1.0000000000000024,
        d=0.5000000000000006,
        e=-0.0,
        f=0.0,
        g=11.430000000000012,
        h=3.0626592694877317,
        j=0.0,
        k=-130.11630454975563,
        name="cy_watercover12",
    )
    surf_45 = openmc.Plane(
        surface_id=45,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-16.3710643889681,
        name="py_watercover18",
    )
    surf_42 = openmc.Plane(
        surface_id=42,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-16.221064388968102,
        name="py_watercover17",
    )
    surf_24 = openmc.Quadric(
        surface_id=24,
        a=0.9330127018922217,
        b=0.06698729810778066,
        c=1.0000000000000024,
        d=0.5000000000000006,
        e=-0.0,
        f=0.0,
        g=11.430000000000012,
        h=3.0626592694877317,
        j=0.0,
        k=-127.55630454975562,
        name="cy_watercover11",
    )
    surf_36 = openmc.Plane(
        surface_id=36,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-3.7210643889681023,
        name="py_watercover15",
    )
    surf_33 = openmc.Plane(
        surface_id=33,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-2.871064388968101,
        name="py_watercover14",
    )
    surf_30 = openmc.Quadric(
        surface_id=30,
        a=0.9330127018922217,
        b=0.06698729810778066,
        c=1.0000000000000024,
        d=0.5000000000000006,
        e=-0.0,
        f=0.0,
        g=11.430000000000012,
        h=3.0626592694877317,
        j=0.0,
        k=-197.5563045497558,
        name="cy_watercover13",
    )
    surf_39 = openmc.Plane(
        surface_id=39,
        a=-0.25881904510252074,
        b=0.9659258262890683,
        c=0.0,
        d=-15.521064388968103,
        name="py_watercover16",
    )
    surf_48 = openmc.Quadric(
        surface_id=48,
        a=0.0669872981077802,
        b=0.9330127018922192,
        c=0.9999999999999993,
        d=-0.49999999999999817,
        e=-0.0,
        f=-0.0,
        g=-5.031993205065669,
        h=18.7796543046465,
        j=0.0,
        k=60.8590928544638,
        name="cox_watercover21",
    )
    surf_51 = openmc.Plane(
        surface_id=51,
        a=0.9659258262890683,
        b=0.25881904510252074,
        c=0.0,
        d=-25.916603371043625,
        name="px_watercover22",
    )
    surf_54 = openmc.Plane(
        surface_id=54,
        a=0.9659258262890683,
        b=0.25881904510252074,
        c=0.0,
        d=-5.916603371043624,
        name="px_watercover23",
    )
    surf_66 = openmc.Quadric(
        surface_id=66,
        a=0.06698729810778066,
        b=0.93301270189222,
        c=1.0000000000000007,
        d=-0.5000000000000001,
        e=-0.0,
        f=-0.0,
        g=0.9726086413127967,
        h=-3.629824865259893,
        j=0.0,
        k=3.227899030403828,
        name="cox_watertube41",
    )
    surf_69 = openmc.Quadric(
        surface_id=69,
        a=0.06698729810778066,
        b=0.93301270189222,
        c=1.0000000000000007,
        d=-0.5000000000000001,
        e=-0.0,
        f=-0.0,
        g=0.9726086413127967,
        h=-3.629824865259893,
        j=0.0,
        k=3.107899030403828,
        name="cox_watertube42",
    )
    surf_72 = openmc.Plane(
        surface_id=72,
        a=0.9659258262890683,
        b=0.25881904510252074,
        c=0.0,
        d=-18.416603371043625,
        name="watertube43",
    )
    surf_75 = openmc.Plane(
        surface_id=75,
        a=0.9659258262890683,
        b=0.25881904510252074,
        c=0.0,
        d=-5.266603371043623,
        name="watertube44",
    )
    surf_80 = openmc.Sphere(
        surface_id=80, x0=0.0, y0=0.0, z0=0.0, r=50.0, name="so_separe61"
    )
    surf_82 = openmc.Plane(
        surface_id=82,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-248.0,
        name="py_rw101",
    )
    surf_84 = openmc.Plane(
        surface_id=84,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=248.0,
        name="py_rw102",
    )
    surf_98 = openmc.Plane(
        surface_id=98, a=1.0, b=0.0, c=0.0, d=-180.0, name="px_rw111"
    )
    surf_100 = openmc.Plane(
        surface_id=100, a=1.0, b=0.0, c=0.0, d=270.0, name="px_rw112"
    )
    surf_114 = openmc.Plane(
        surface_id=114,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=-248.0,
        name="pz_rw121",
    )
    surf_120 = openmc.Plane(
        surface_id=120,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=250.0,
        name="pz_rw124",
    )
    surf_86 = openmc.Plane(
        surface_id=86,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-250.0,
        name="py_rw103",
    )
    surf_88 = openmc.Plane(
        surface_id=88,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=250.0,
        name="py_rw104",
    )
    surf_102 = openmc.Plane(
        surface_id=102, a=1.0, b=0.0, c=0.0, d=-182.0, name="px_rw113"
    )
    surf_104 = openmc.Plane(
        surface_id=104, a=1.0, b=0.0, c=0.0, d=272.0, name="px_rw114"
    )
    surf_118 = openmc.Plane(
        surface_id=118,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=-250.0,
        name="pz_rw123",
    )
    surf_124 = openmc.Plane(
        surface_id=124,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=259.0,
        name="pz_rw126",
    )
    surf_90 = openmc.Plane(
        surface_id=90,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-259.0,
        name="py_rw105",
    )
    surf_92 = openmc.Plane(
        surface_id=92,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=259.0,
        name="py_rw106",
    )
    surf_106 = openmc.Plane(
        surface_id=106, a=1.0, b=0.0, c=0.0, d=-191.0, name="px_rw115"
    )
    surf_108 = openmc.Plane(
        surface_id=108, a=1.0, b=0.0, c=0.0, d=281.0, name="px_rw116"
    )
    surf_122 = openmc.Plane(
        surface_id=122,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=-259.0,
        name="pz_rw125",
    )
    surf_128 = openmc.Plane(
        surface_id=128,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=268.0,
        name="pz_rw128",
    )
    surf_94 = openmc.Plane(
        surface_id=94,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-268.0,
        boundary_type="vacuum",
        name="py_rw107",
    )
    surf_96 = openmc.Plane(
        surface_id=96,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=268.0,
        boundary_type="vacuum",
        name="py_rw108",
    )
    surf_110 = openmc.Plane(
        surface_id=110,
        a=1.0,
        b=0.0,
        c=0.0,
        d=-200.0,
        boundary_type="vacuum",
        name="px_rw117",
    )
    surf_112 = openmc.Plane(
        surface_id=112,
        a=1.0,
        b=0.0,
        c=0.0,
        d=290.0,
        boundary_type="vacuum",
        name="px_rw118",
    )
    surf_126 = openmc.Plane(
        surface_id=126,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=-268.0,
        boundary_type="vacuum",
        name="pz_rw127",
    )
    surf_134 = openmc.Plane(
        surface_id=134,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-126.0,
        name="py_osa133",
    )
    surf_136 = openmc.Plane(
        surface_id=136,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=126.0,
        name="py_osa134",
    )
    surf_150 = openmc.Plane(
        surface_id=150, a=1.0, b=0.0, c=0.0, d=131.2, name="px_osa143"
    )
    surf_116 = openmc.Plane(
        surface_id=116,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=248.0,
        name="pz_rw122",
    )
    surf_138 = openmc.Plane(
        surface_id=138,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-134.2,
        name="py_osa135",
    )
    surf_140 = openmc.Plane(
        surface_id=140,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=134.2,
        name="py_osa136",
    )
    surf_152 = openmc.Plane(
        surface_id=152, a=1.0, b=0.0, c=0.0, d=139.2, name="px_osa144"
    )
    surf_130 = openmc.Plane(
        surface_id=130,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-125.0,
        name="py_osa131",
    )
    surf_132 = openmc.Plane(
        surface_id=132,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=125.0,
        name="py_osa132",
    )
    surf_148 = openmc.Plane(
        surface_id=148, a=1.0, b=0.0, c=0.0, d=130.0, name="px_osa142"
    )
    surf_158 = openmc.Plane(
        surface_id=158,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=450.0,
        name="pz_osa152",
    )
    surf_142 = openmc.Plane(
        surface_id=142,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-142.2,
        name="py_osa137",
    )
    surf_144 = openmc.Plane(
        surface_id=144,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=142.2,
        name="py_osa138",
    )
    surf_154 = openmc.Plane(
        surface_id=154, a=1.0, b=0.0, c=0.0, d=147.2, name="px_osa145"
    )
    surf_146 = openmc.Plane(
        surface_id=146, a=1.0, b=0.0, c=0.0, d=-176.5, name="px_osa141"
    )
    surf_156 = openmc.Plane(
        surface_id=156,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=120.0,
        name="pz_osa151",
    )
    surf_161 = openmc.YPlane(surface_id=161, y0=-50.0, name="py_osa160")
    surf_163 = openmc.XPlane(surface_id=163, x0=30.0, name="px_osa162")
    surf_164 = openmc.ZPlane(surface_id=164, z0=-30.0, name="pz_osa163")
    surf_165 = openmc.ZPlane(surface_id=165, z0=30.0, name="pz_osa164")
    surf_162 = openmc.YPlane(surface_id=162, y0=100.0, name="py_osa161")
    surf_225 = openmc.Plane(
        surface_id=225,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=149.9996537907726,
        name="pz_sas231",
    )
    surf_264 = openmc.Plane(
        surface_id=264,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-2.0001320915440957,
        name="py_sas271",
    )
    surf_267 = openmc.Plane(
        surface_id=267,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=1.9998679084559043,
        name="py_sas272",
    )
    surf_282 = openmc.Plane(
        surface_id=282, a=1.0, b=0.0, c=0.0, d=-2.0, name="px_sas285"
    )
    surf_285 = openmc.Plane(
        surface_id=285, a=1.0, b=0.0, c=0.0, d=2.0, name="px_sas286"
    )
    surf_300 = openmc.Plane(
        surface_id=300,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=150.9996537907726,
        name="pz_sas291",
    )
    surf_174 = openmc.Plane(
        surface_id=174,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-50.000132091544096,
        name="py_sas207",
    )
    surf_177 = openmc.Plane(
        surface_id=177,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=54.999867908455904,
        name="py_sas208",
    )
    surf_201 = openmc.Plane(
        surface_id=201, a=1.0, b=0.0, c=0.0, d=-70.0, name="px_sas224"
    )
    surf_207 = openmc.Plane(
        surface_id=207, a=1.0, b=0.0, c=0.0, d=55.0, name="px_sas227"
    )
    surf_228 = openmc.Plane(
        surface_id=228,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=169.9996537907726,
        name="pz_sas235",
    )
    surf_192 = openmc.Plane(
        surface_id=192, a=1.0, b=0.0, c=0.0, d=-146.5, name="px_sas221"
    )
    surf_234 = openmc.Plane(
        surface_id=234,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=199.9996537907726,
        name="pz_sas238",
    )
    surf_243 = openmc.Plane(
        surface_id=243,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=234.9996537907726,
        name="pz_sas243",
    )
    surf_249 = openmc.Plane(
        surface_id=249,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=299.99965379077264,
        name="pz_sas245",
    )
    surf_252 = openmc.Plane(
        surface_id=252,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=349.99965379077264,
        name="pz_sas250",
    )
    surf_186 = openmc.Plane(
        surface_id=186,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-75.00013209154409,
        name="py_sas211",
    )
    surf_189 = openmc.Plane(
        surface_id=189,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=74.99986790845591,
        name="py_sas212",
    )
    surf_195 = openmc.Plane(
        surface_id=195, a=1.0, b=0.0, c=0.0, d=-90.0, name="px_sas222"
    )
    surf_210 = openmc.Plane(
        surface_id=210, a=1.0, b=0.0, c=0.0, d=75.0, name="px_sas228"
    )
    surf_180 = openmc.Plane(
        surface_id=180,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-70.00013209154409,
        name="py_sas209",
    )
    surf_183 = openmc.Plane(
        surface_id=183,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=69.99986790845591,
        name="py_sas210",
    )
    surf_198 = openmc.Plane(
        surface_id=198, a=1.0, b=0.0, c=0.0, d=-85.0, name="px_sas223"
    )
    surf_219 = openmc.Plane(
        surface_id=219, a=1.0, b=0.0, c=0.0, d=-15.0, name="px_sas263"
    )
    surf_222 = openmc.Plane(
        surface_id=222, a=1.0, b=0.0, c=0.0, d=15.0, name="px_sas264"
    )
    surf_168 = openmc.Plane(
        surface_id=168,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=-15.000132091544096,
        name="py_sas203",
    )
    surf_171 = openmc.Plane(
        surface_id=171,
        a=0.0,
        b=0.984807753012208,
        c=-0.17364817766693033,
        d=14.999867908455904,
        name="py_sas204",
    )
    surf_261 = openmc.Plane(
        surface_id=261,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=269.99965379077264,
        name="pz_sas253",
    )
    surf_216 = openmc.Plane(
        surface_id=216, a=1.0, b=0.0, c=0.0, d=-45.0, name="px_sas262"
    )
    surf_213 = openmc.Plane(
        surface_id=213, a=1.0, b=0.0, c=0.0, d=-75.0, name="px_sas261"
    )
    surf_270 = openmc.Plane(
        surface_id=270, a=1.0, b=0.0, c=0.0, d=-62.0, name="px_sas281"
    )
    surf_273 = openmc.Plane(
        surface_id=273, a=1.0, b=0.0, c=0.0, d=-58.0, name="px_sas282"
    )
    surf_324 = openmc.Plane(
        surface_id=324,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=350.99965379077264,
        name="pz_sas299",
    )
    surf_276 = openmc.Plane(
        surface_id=276, a=1.0, b=0.0, c=0.0, d=-32.0, name="px_sas283"
    )
    surf_279 = openmc.Plane(
        surface_id=279, a=1.0, b=0.0, c=0.0, d=-28.0, name="px_sas284"
    )
    surf_288 = openmc.Plane(
        surface_id=288, a=1.0, b=0.0, c=0.0, d=28.0, name="px_sas287"
    )
    surf_291 = openmc.Plane(
        surface_id=291, a=1.0, b=0.0, c=0.0, d=32.0, name="px_sas288"
    )
    surf_294 = openmc.Plane(
        surface_id=294, a=1.0, b=0.0, c=0.0, d=58.0, name="px_sas289"
    )
    surf_297 = openmc.Plane(
        surface_id=297, a=1.0, b=0.0, c=0.0, d=62.0, name="px_sas290"
    )
    surf_303 = openmc.Plane(
        surface_id=303,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=170.9996537907726,
        name="pz_sas292",
    )
    surf_306 = openmc.Plane(
        surface_id=306,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=219.4996537907726,
        name="pz_sas293",
    )
    surf_309 = openmc.Plane(
        surface_id=309,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=220.4996537907726,
        name="pz_sas294",
    )
    surf_312 = openmc.Plane(
        surface_id=312,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=284.49965379077264,
        name="pz_sas295",
    )
    surf_315 = openmc.Plane(
        surface_id=315,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=285.49965379077264,
        name="pz_sas296",
    )
    surf_318 = openmc.Plane(
        surface_id=318,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=309.49965379077264,
        name="pz_sas297",
    )
    surf_321 = openmc.Plane(
        surface_id=321,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=310.49965379077264,
        name="pz_sas298",
    )
    surf_160 = openmc.Plane(
        surface_id=160,
        a=0.0,
        b=0.17364817766693033,
        c=0.984807753012208,
        d=500.0,
        boundary_type="vacuum",
        name="pz_osa155",
    )

    # regions
    region_1 = +surf_3 & -surf_6 & -surf_18 & +surf_21
    region_2 = +surf_6 & -surf_9 & -surf_18 & +surf_21
    region_3 = +surf_76 & -surf_77 & -surf_78 & +surf_79
    region_4 = (
        -surf_3 & -surf_15 & +
        surf_57 & ~(+surf_76 & -surf_77 & -surf_78 & +surf_79)
    )
    region_5 = +surf_15 & -surf_18 & +surf_57 & -surf_3
    region_6 = -surf_63 & +surf_60 & -surf_57
    region_11 = +surf_9 & -surf_12 & -surf_27 & +surf_21
    region_12 = +surf_45 & -surf_42 & +surf_18 & -surf_27
    region_13 = +surf_24 & -surf_27 & +surf_42 & -surf_9
    region_14 = +surf_36 & -surf_33 & +surf_27 & -surf_30
    region_15 = +surf_45 & -surf_39 & +surf_27 & -surf_30
    region_16 = +surf_18 & -surf_24 & +surf_42 & -surf_9
    region_17 = -surf_48 & +surf_51 & +surf_27 & -surf_54
    region_21 = +surf_66 & -surf_69 & +surf_72 & -surf_75
    region_22 = -surf_66 & +surf_72 & -surf_75
    region_31 = (
        (-surf_80 & -surf_60)
        | (-surf_80 & +surf_60 & -surf_57 & +surf_63)
        | (-surf_80 & +surf_57 & -surf_45 & +surf_18)
    )
    region_32 = -surf_80 & +surf_45 & -surf_39 & +surf_30
    region_33 = (
        ~(-surf_48 & +surf_51 & +surf_27 & -surf_54)
        & -surf_80
        & +surf_39
        & -surf_36
        & +surf_27
    )
    region_34 = -surf_80 & +surf_36 & -surf_33 & +surf_30
    region_35 = (
        -surf_80 & -surf_27 & +
        surf_12 & ~(-surf_69 & +surf_72 & -surf_75) & +surf_21
    )
    region_36 = -surf_80 & +surf_27 & +surf_33
    region_101 = (
        ~(+surf_82 & -surf_84 & +surf_98 & -surf_100 & +surf_114 & -surf_120)
        & +surf_86
        & -surf_88
        & +surf_102
        & -surf_104
        & +surf_118
        & -surf_120
    )
    region_102 = (
        ~(+surf_86 & -surf_88 & +surf_102 & -surf_104 & +surf_118 & -surf_124)
        & +surf_90
        & -surf_92
        & +surf_106
        & -surf_108
        & +surf_122
        & -surf_124
    )
    region_103 = (
        ~(+surf_90 & -surf_92 & +surf_106 & -surf_108 & +surf_122 & -surf_128)
        & +surf_94
        & -surf_96
        & +surf_110
        & -surf_112
        & +surf_126
        & -surf_128
    )
    region_111 = (
        ~(+surf_134 & -surf_136 & +surf_98 & -surf_150 & +surf_116 & -surf_120)
        & +surf_82
        & -surf_84
        & +surf_98
        & -surf_100
        & +surf_116
        & -surf_120
    )
    region_112 = (
        ~(+surf_134 & -surf_136 & +surf_98 & -surf_150 & +surf_120 & -surf_124)
        & +surf_86
        & -surf_88
        & +surf_102
        & -surf_104
        & +surf_120
        & -surf_124
    )
    region_113 = (
        ~(+surf_138 & -surf_140 & +surf_106 & -surf_152 & +surf_124 & -surf_128)
        & +surf_90
        & -surf_92
        & +surf_106
        & -surf_108
        & +surf_124
        & -surf_128
    )
    region_121 = (
        ~(+surf_130 & -surf_132 & +surf_98 & -surf_148 & +surf_116 & -surf_158)
        & +surf_134
        & -surf_136
        & +surf_98
        & -surf_150
        & +surf_116
        & -surf_158
    )
    region_122 = (
        ~(+surf_134 & -surf_136 & +surf_98 & -surf_150 & +surf_124 & -surf_158)
        & +surf_138
        & -surf_140
        & +surf_106
        & -surf_152
        & +surf_124
        & -surf_158
    )
    region_123 = (
        ~(+surf_138 & -surf_140 & +surf_106 & -surf_152 & +surf_128 & -surf_158)
        & +surf_142
        & -surf_144
        & +surf_110
        & -surf_154
        & +surf_128
        & -surf_158
    )
    region_131 = +surf_98 & -surf_146 & +surf_130 & -surf_132 & +surf_156 & -surf_158
    region_132 = +surf_82 & -surf_161 & +surf_98 & -surf_163 & +surf_164 & -surf_165
    region_133 = +surf_162 & -surf_84 & +surf_98 & -surf_163 & +surf_164 & -surf_165
    region_135 = (
        ~(+surf_82 & -surf_161 & +surf_98 & -surf_163 & +surf_164 & -surf_165)
        & ~(+surf_162 & -surf_84 & +surf_98 & -surf_163 & +surf_164 & -surf_165)
        & +surf_82
        & -surf_84
        & +surf_98
        & -surf_100
        & +surf_114
        & -surf_156
        & +surf_80
    )
    region_136 = (
        ~(+surf_130 & -surf_132 & +surf_98 & -surf_148 & +surf_156 & -surf_116)
        & +surf_82
        & -surf_84
        & +surf_98
        & -surf_100
        & +surf_156
        & -surf_116
    )
    region_201 = +surf_130 & -surf_132 & + \
        surf_146 & -surf_148 & +surf_156 & -surf_225
    region_202 = (
        ~(+surf_264 & -surf_267 & +surf_282 & -surf_285 & +surf_225 & -surf_300)
        & +surf_174
        & -surf_177
        & +surf_201
        & -surf_207
        & +surf_225
        & -surf_228
    )
    region_211 = +surf_130 & -surf_132 & + \
        surf_146 & -surf_192 & +surf_225 & -surf_234
    region_212 = +surf_130 & -surf_132 & + \
        surf_146 & -surf_192 & +surf_234 & -surf_243
    region_213 = +surf_130 & -surf_132 & + \
        surf_146 & -surf_192 & +surf_243 & -surf_249
    region_214 = +surf_130 & -surf_132 & + \
        surf_146 & -surf_192 & +surf_249 & -surf_252
    region_221 = (
        ~(+surf_174 & -surf_177 & +surf_201 & -surf_207 & +surf_225 & -surf_228)
        & +surf_130
        & -surf_132
        & +surf_192
        & -surf_148
        & +surf_225
        & -surf_228
    )
    region_222 = (
        ~(+surf_186 & -surf_189 & +surf_195 & -surf_210 & +surf_228 & -surf_234)
        & +surf_130
        & -surf_132
        & +surf_192
        & -surf_148
        & +surf_228
        & -surf_234
    )
    region_223 = (
        ~(+surf_186 & -surf_189 & +surf_195 & -surf_210 & +surf_234 & -surf_243)
        & +surf_130
        & -surf_132
        & +surf_192
        & -surf_148
        & +surf_234
        & -surf_243
    )
    region_224 = (
        ~(+surf_186 & -surf_189 & +surf_195 & -surf_210 & +surf_243 & -surf_249)
        & +surf_130
        & -surf_132
        & +surf_192
        & -surf_148
        & +surf_243
        & -surf_249
    )
    region_225 = (
        ~(+surf_186 & -surf_189 & +surf_195 & -surf_210 & +surf_249 & -surf_252)
        & +surf_130
        & -surf_132
        & +surf_192
        & -surf_148
        & +surf_249
        & -surf_252
    )
    region_231 = (
        ~(+surf_180 & -surf_183 & +surf_198 & -surf_210 & +surf_228 & -surf_234)
        & +surf_186
        & -surf_189
        & +surf_195
        & -surf_210
        & +surf_228
        & -surf_234
    )
    region_232 = (
        ~(+surf_180 & -surf_183 & +surf_198 & -surf_210 & +surf_234 & -surf_243)
        & +surf_186
        & -surf_189
        & +surf_195
        & -surf_210
        & +surf_234
        & -surf_243
    )
    region_233 = (
        ~(+surf_180 & -surf_183 & +surf_198 & -surf_210 & +surf_243 & -surf_252)
        & +surf_186
        & -surf_189
        & +surf_195
        & -surf_210
        & +surf_243
        & -surf_252
    )
    region_241 = (
        ~(+surf_219 & -surf_222 & +surf_168 & -surf_171 & +surf_228 & -surf_234)
        & +surf_198
        & -surf_210
        & +surf_180
        & -surf_183
        & +surf_228
        & -surf_234
    )
    region_242 = (
        ~(+surf_219 & -surf_222 & +surf_168 & -surf_171 & +surf_234 & -surf_243)
        & +surf_198
        & -surf_210
        & +surf_180
        & -surf_183
        & +surf_234
        & -surf_243
    )
    region_243 = (
        ~(+surf_219 & -surf_222 & +surf_168 & -surf_171 & +surf_243 & -surf_261)
        & ~(+surf_219 & -surf_222 & +surf_168 & -surf_171 & +surf_261 & -surf_249)
        & ~(+surf_216 & -surf_219 & +surf_168 & -surf_171 & +surf_261 & -surf_249)
        & ~(+surf_213 & -surf_216 & +surf_168 & -surf_171 & +surf_261 & -surf_249)
        & +surf_198
        & -surf_210
        & +surf_180
        & -surf_183
        & +surf_243
        & -surf_249
    )
    region_244 = (
        ~(+surf_213 & -surf_216 & +surf_168 & -surf_171 & +surf_249 & -surf_252)
        & +surf_198
        & -surf_210
        & +surf_180
        & -surf_183
        & +surf_249
        & -surf_252
    )
    region_253 = (
        ~(+surf_270 & -surf_273 & +surf_264 & -surf_267 & +surf_252 & -surf_324)
        & ~(+surf_276 & -surf_279 & +surf_264 & -surf_267 & +surf_252 & -surf_324)
        & ~(+surf_282 & -surf_285 & +surf_264 & -surf_267 & +surf_252 & -surf_324)
        & ~(+surf_288 & -surf_291 & +surf_264 & -surf_267 & +surf_252 & -surf_324)
        & ~(+surf_294 & -surf_297 & +surf_264 & -surf_267 & +surf_252 & -surf_324)
        & +surf_146
        & -surf_148
        & +surf_130
        & -surf_132
        & +surf_252
        & -surf_158
    )
    region_300 = +surf_264 & -surf_267 & + \
        surf_282 & -surf_285 & +surf_225 & -surf_300
    region_301 = +surf_264 & -surf_267 & + \
        surf_282 & -surf_285 & +surf_228 & -surf_303
    region_302 = +surf_264 & -surf_267 & + \
        surf_282 & -surf_285 & +surf_306 & -surf_309
    region_303 = +surf_264 & -surf_267 & + \
        surf_282 & -surf_285 & +surf_312 & -surf_315
    region_304 = +surf_264 & -surf_267 & + \
        surf_276 & -surf_279 & +surf_312 & -surf_315
    region_305 = +surf_264 & -surf_267 & + \
        surf_270 & -surf_273 & +surf_312 & -surf_315
    region_306 = +surf_264 & -surf_267 & + \
        surf_270 & -surf_273 & +surf_318 & -surf_321
    region_307 = +surf_264 & -surf_267 & + \
        surf_270 & -surf_273 & +surf_252 & -surf_324
    region_308 = +surf_264 & -surf_267 & + \
        surf_276 & -surf_279 & +surf_252 & -surf_324
    region_309 = +surf_264 & -surf_267 & + \
        surf_282 & -surf_285 & +surf_252 & -surf_324
    region_310 = +surf_264 & -surf_267 & + \
        surf_288 & -surf_291 & +surf_252 & -surf_324
    region_311 = +surf_264 & -surf_267 & + \
        surf_294 & -surf_297 & +surf_252 & -surf_324
    region_401 = (
        ~(+surf_264 & -surf_267 & +surf_282 & -surf_285 & +surf_228 & -surf_303)
        & +surf_168
        & -surf_171
        & +surf_219
        & -surf_222
        & +surf_228
        & -surf_234
    )
    region_402 = (
        ~(+surf_264 & -surf_267 & +surf_282 & -surf_285 & +surf_306 & -surf_309)
        & +surf_168
        & -surf_171
        & +surf_219
        & -surf_222
        & +surf_234
        & -surf_243
    )
    region_403 = +surf_168 & -surf_171 & + \
        surf_219 & -surf_222 & +surf_243 & -surf_261
    region_404 = (
        ~(+surf_264 & -surf_267 & +surf_282 & -surf_285 & +surf_312 & -surf_315)
        & +surf_168
        & -surf_171
        & +surf_219
        & -surf_222
        & +surf_261
        & -surf_249
    )
    region_405 = (
        ~(+surf_264 & -surf_267 & +surf_276 & -surf_279 & +surf_312 & -surf_315)
        & +surf_168
        & -surf_171
        & +surf_216
        & -surf_219
        & +surf_261
        & -surf_249
    )
    region_406 = (
        ~(+surf_264 & -surf_267 & +surf_270 & -surf_273 & +surf_312 & -surf_315)
        & +surf_168
        & -surf_171
        & +surf_213
        & -surf_216
        & +surf_261
        & -surf_249
    )
    region_407 = (
        ~(+surf_264 & -surf_267 & +surf_270 & -surf_273 & +surf_318 & -surf_321)
        & +surf_168
        & -surf_171
        & +surf_213
        & -surf_216
        & +surf_249
        & -surf_252
    )
    region_502 = (
        +surf_128
        & -surf_160
        & ~(+surf_142 & -surf_144 & +surf_110 & -surf_154 & +surf_128 & -surf_158)
        & ~(+surf_96 | -surf_94 | +surf_112 | -surf_110 | -surf_126)
    )

    # cells
    cell_1 = openmc.Cell(cell_id=1, region=region_1, fill=mat_1)
    cell_2 = openmc.Cell(cell_id=2, region=region_2, fill=mat_2)
    cell_3 = openmc.Cell(cell_id=3, region=region_3, fill=mat_3)
    cell_4 = openmc.Cell(cell_id=4, region=region_4, fill=None)
    cell_5 = openmc.Cell(cell_id=5, region=region_5, fill=mat_4)
    cell_6 = openmc.Cell(cell_id=6, region=region_6, fill=mat_5)
    cell_11 = openmc.Cell(cell_id=11, region=region_11, fill=mat_4)
    cell_12 = openmc.Cell(cell_id=12, region=region_12, fill=mat_4)
    cell_13 = openmc.Cell(cell_id=13, region=region_13, fill=mat_5)
    cell_14 = openmc.Cell(cell_id=14, region=region_14, fill=mat_5)
    cell_15 = openmc.Cell(cell_id=15, region=region_15, fill=mat_5)
    cell_16 = openmc.Cell(cell_id=16, region=region_16, fill=mat_21)
    cell_17 = openmc.Cell(cell_id=17, region=region_17, fill=mat_22)
    cell_21 = openmc.Cell(cell_id=21, region=region_21, fill=mat_4)
    cell_22 = openmc.Cell(cell_id=22, region=region_22, fill=mat_2)
    cell_31 = openmc.Cell(cell_id=31, region=region_31, fill=mat_6)
    cell_32 = openmc.Cell(cell_id=32, region=region_32, fill=mat_6)
    cell_33 = openmc.Cell(cell_id=33, region=region_33, fill=mat_6)
    cell_34 = openmc.Cell(cell_id=34, region=region_34, fill=mat_6)
    cell_35 = openmc.Cell(cell_id=35, region=region_35, fill=mat_6)
    cell_36 = openmc.Cell(cell_id=36, region=region_36, fill=mat_6)
    cell_101 = openmc.Cell(cell_id=101, region=region_101, fill=mat_7)
    cell_102 = openmc.Cell(cell_id=102, region=region_102, fill=mat_8)
    cell_103 = openmc.Cell(cell_id=103, region=region_103, fill=mat_8)
    cell_111 = openmc.Cell(cell_id=111, region=region_111, fill=mat_7)
    cell_112 = openmc.Cell(cell_id=112, region=region_112, fill=mat_8)
    cell_113 = openmc.Cell(cell_id=113, region=region_113, fill=mat_8)
    cell_121 = openmc.Cell(cell_id=121, region=region_121, fill=mat_11)
    cell_122 = openmc.Cell(cell_id=122, region=region_122, fill=mat_8)
    cell_123 = openmc.Cell(cell_id=123, region=region_123, fill=mat_8)
    cell_131 = openmc.Cell(cell_id=131, region=region_131, fill=mat_11)
    cell_132 = openmc.Cell(cell_id=132, region=region_132, fill=mat_41)
    cell_133 = openmc.Cell(cell_id=133, region=region_133, fill=mat_42)
    cell_135 = openmc.Cell(cell_id=135, region=region_135, fill=mat_6)
    cell_136 = openmc.Cell(cell_id=136, region=region_136, fill=mat_6)
    cell_201 = openmc.Cell(cell_id=201, region=region_201, fill=mat_6)
    cell_202 = openmc.Cell(cell_id=202, region=region_202, fill=mat_6)
    cell_211 = openmc.Cell(cell_id=211, region=region_211, fill=mat_11)
    cell_212 = openmc.Cell(cell_id=212, region=region_212, fill=mat_11)
    cell_213 = openmc.Cell(cell_id=213, region=region_213, fill=mat_11)
    cell_214 = openmc.Cell(cell_id=214, region=region_214, fill=mat_11)
    cell_221 = openmc.Cell(cell_id=221, region=region_221, fill=mat_12)
    cell_222 = openmc.Cell(cell_id=222, region=region_222, fill=mat_12)
    cell_223 = openmc.Cell(cell_id=223, region=region_223, fill=mat_12)
    cell_224 = openmc.Cell(cell_id=224, region=region_224, fill=mat_12)
    cell_225 = openmc.Cell(cell_id=225, region=region_225, fill=mat_12)
    cell_231 = openmc.Cell(cell_id=231, region=region_231, fill=mat_11)
    cell_232 = openmc.Cell(cell_id=232, region=region_232, fill=mat_11)
    cell_233 = openmc.Cell(cell_id=233, region=region_233, fill=mat_11)
    cell_241 = openmc.Cell(cell_id=241, region=region_241, fill=mat_11)
    cell_242 = openmc.Cell(cell_id=242, region=region_242, fill=mat_11)
    cell_243 = openmc.Cell(cell_id=243, region=region_243, fill=mat_11)
    cell_244 = openmc.Cell(cell_id=244, region=region_244, fill=mat_11)
    cell_253 = openmc.Cell(cell_id=253, region=region_253, fill=mat_6)
    cell_300 = openmc.Cell(cell_id=300, region=region_300, fill=mat_6)
    cell_301 = openmc.Cell(
        cell_id=301, region=region_301, fill=mat_9)  # detector
    cell_302 = openmc.Cell(
        cell_id=302, region=region_302, fill=mat_9)  # detector
    cell_303 = openmc.Cell(
        cell_id=303, region=region_303, fill=mat_9)  # detector
    cell_304 = openmc.Cell(
        cell_id=304, region=region_304, fill=mat_9)  # detector
    cell_305 = openmc.Cell(
        cell_id=305, region=region_305, fill=mat_9)  # detector
    cell_306 = openmc.Cell(
        cell_id=306, region=region_306, fill=mat_9)  # detector
    cell_307 = openmc.Cell(
        cell_id=307, region=region_307, fill=mat_9)  # detector
    cell_308 = openmc.Cell(
        cell_id=308, region=region_308, fill=mat_9)  # detector
    cell_309 = openmc.Cell(
        cell_id=309, region=region_309, fill=mat_9)  # detector
    cell_310 = openmc.Cell(
        cell_id=310, region=region_310, fill=mat_9)  # detector
    cell_311 = openmc.Cell(
        cell_id=311, region=region_311, fill=mat_9)  # detector
    cell_401 = openmc.Cell(cell_id=401, region=region_401, fill=mat_6)
    cell_402 = openmc.Cell(cell_id=402, region=region_402, fill=mat_6)
    cell_403 = openmc.Cell(cell_id=403, region=region_403, fill=mat_6)
    cell_404 = openmc.Cell(cell_id=404, region=region_404, fill=mat_6)
    cell_405 = openmc.Cell(cell_id=405, region=region_405, fill=mat_6)
    cell_406 = openmc.Cell(cell_id=406, region=region_406, fill=mat_6)
    cell_407 = openmc.Cell(cell_id=407, region=region_407, fill=mat_6)
    cell_502 = openmc.Cell(cell_id=502, region=region_502, fill=None)

    # create root universe
    universe = openmc.Universe(
        cells=[
            cell_1,
            cell_2,
            cell_3,
            cell_4,
            cell_5,
            cell_6,
            cell_11,
            cell_12,
            cell_13,
            cell_14,
            cell_15,
            cell_16,
            cell_17,
            cell_21,
            cell_22,
            cell_31,
            cell_32,
            cell_33,
            cell_34,
            cell_35,
            cell_36,
            cell_101,
            cell_102,
            cell_103,
            cell_111,
            cell_112,
            cell_113,
            cell_121,
            cell_122,
            cell_123,
            cell_131,
            cell_132,
            cell_133,
            cell_135,
            cell_136,
            cell_201,
            cell_202,
            cell_211,
            cell_212,
            cell_213,
            cell_214,
            cell_221,
            cell_222,
            cell_223,
            cell_224,
            cell_225,
            cell_231,
            cell_232,
            cell_233,
            cell_241,
            cell_242,
            cell_243,
            cell_244,
            cell_253,
            cell_300,
            cell_301,
            cell_302,
            cell_303,
            cell_304,
            cell_305,
            cell_306,
            cell_307,
            cell_308,
            cell_309,
            cell_310,
            cell_311,
            cell_401,
            cell_402,
            cell_403,
            cell_404,
            cell_405,
            cell_406,
            cell_407,
            cell_502,
        ]
    )

    # create geometry instance
    model.geometry = openmc.Geometry(universe)

    ###############################################################################
    # Define problem settings

    # weight windows from wwinps
    ww = openmc.wwinp_to_wws("weight_windows.cadis.wwinp")

    # Indicate how many particles to run
    settings = openmc.Settings(run_mode='fixed source')
    settings.batches = args.batches
    settings.particles = args.particles
    # fng source
    fng2fns_source = fng_source(
        center=(0, 0, 0), reference_uvw=(0, 1, 0), beam_energy=230)
    settings.source = fng2fns_source
    settings.weight_windows = ww
    settings.output = {'tallies': False}

    ############################################################################
    # define tallies

    model.tallies = openmc.Tallies()

    # filters
    # particle filters
    neutron_filter = openmc.ParticleFilter("neutron")
    particle_filter = openmc.ParticleFilter(["neutron", "photon"])
    # cell filters
    # detector00 is never used and it is not part of the experiment
    detector_cell_filter = openmc.CellFilter(
        [
            cell_301,
            cell_302,
            cell_303,
            cell_304,
            cell_305,
            cell_306,
            cell_307,
            cell_308,
            cell_309,
            cell_310,
            cell_311,
        ]
    )

    # energy filters
    # Neutron Spectrum in 125-Energy Bin
    neutron_energy_filter = openmc.EnergyFilter(
        1e6*np.array(
            [
                1.0010e-11,
                3.2241e-07,
                5.3156e-07,
                8.7640e-07,
                1.4449e-06,
                2.3823e-06,
                3.9278e-06,
                6.4758e-06,
                1.0677e-05,
                1.7603e-05,
                2.9023e-05,
                4.7850e-05,
                7.8891e-05,
                1.3007e-04,
                2.1445e-04,
                3.5357e-04,
                5.8293e-04,
                9.6110e-04,
                1.2341e-03,
                1.5846e-03,
                2.0346e-03,
                2.6125e-03,
                3.3546e-03,
                4.3073e-03,
                5.5307e-03,
                7.1016e-03,
                9.1186e-03,
                1.1709e-02,
                1.5034e-02,
                1.9304e-02,
                2.1874e-02,
                2.4787e-02,
                2.8087e-02,
                3.1827e-02,
                3.6065e-02,
                4.0867e-02,
                4.6308e-02,
                5.2474e-02,
                5.9461e-02,
                6.7378e-02,
                7.6349e-02,
                8.6515e-02,
                9.8035e-02,
                1.1109e-01,
                1.2588e-01,
                1.4264e-01,
                1.6163e-01,
                1.8315e-01,
                2.0754e-01,
                2.3517e-01,
                2.6649e-01,
                3.0197e-01,
                3.4217e-01,
                3.8774e-01,
                4.3936e-01,
                4.9786e-01,
                5.6415e-01,
                6.3927e-01,
                7.2438e-01,
                8.2084e-01,
                9.3013e-01,
                1.0540e00,
                1.1943e00,
                1.3533e00,
                1.5335e00,
                1.7377e00,
                1.8498e00,
                1.9691e00,
                2.0961e00,
                2.2313e00,
                2.3752e00,
                2.5284e00,
                2.6914e00,
                2.8650e00,
                3.0498e00,
                3.2465e00,
                3.4559e00,
                3.6787e00,
                3.9160e00,
                4.1686e00,
                4.4374e00,
                4.7236e00,
                5.0282e00,
                5.3525e00,
                5.6978e00,
                6.0652e00,
                6.4564e00,
                6.8728e00,
                7.3161e00,
                7.7879e00,
                8.2902e00,
                8.8249e00,
                9.3940e00,
                9.9999e00,
                1.0157e01,
                1.0317e01,
                1.0480e01,
                1.0645e01,
                1.0812e01,
                1.0983e01,
                1.1156e01,
                1.1331e01,
                1.1510e01,
                1.1691e01,
                1.1875e01,
                1.2062e01,
                1.2252e01,
                1.2445e01,
                1.2641e01,
                1.2840e01,
                1.3042e01,
                1.3248e01,
                1.3456e01,
                1.3668e01,
                1.3883e01,
                1.4102e01,
                1.4324e01,
                1.4550e01,
                1.4779e01,
                1.5012e01,
                1.5248e01,
                1.5488e01,
                1.5732e01,
            ]
        )
    )  # mcnp uses MeV, openmc uses eV

    # dosimetry tallies from IRDFF-II nuclear data library
    nb93_n2n_acef = irdff.path + "dos-irdff2-4125.acef"
    in115_nn_acef = irdff.path + "dos-irdff2-4931.acef"
    au197_ng_acef = irdff.path + "dos-irdff2-7925_modified.acef"
    irdff_xs = [nb93_n2n_acef, in115_nn_acef, au197_ng_acef]
    reactions = [11016, 11004, 102]
    nuclides = ['nb93', 'in115', 'au197']

    for n, r, x in zip(nuclides, reactions, irdff_xs):
        # onaxis1 tally
        tally1 = openmc.Tally(name=f"rr_{n}")
        irdff_xs = irdff.cross_section(x)
        multiplier = openmc.EnergyFunctionFilter.from_tabulated1d(
            irdff_xs[r])
        tally1.filters = [detector_cell_filter, neutron_filter, multiplier]
        tally1.scores = ["flux"]
        model.tallies.extend([tally1])

    # neutrons energy spectrum tallies - same energy bins as mcnp model
    # cell tally - neutron spectrum at detector
    det_cell = [cell_303, cell_305, cell_307, cell_309]
    det_pos = ['3', '5', '7', '9']
    for dc, dp in zip(det_cell, det_pos):
        cell_filter = openmc.CellFilter([dc])
        tally101 = openmc.Tally(name=f"nspectrum_detector{dp}")
        tally101.filters = [cell_filter, neutron_filter, neutron_energy_filter]
        tally101.scores = ["flux"]
        model.tallies.extend([tally101])

    cwd = 'results'
    model.settings = settings

    return model.run(cwd=cwd, threads=args.threads)


if __name__ == "__main__":
    main()
