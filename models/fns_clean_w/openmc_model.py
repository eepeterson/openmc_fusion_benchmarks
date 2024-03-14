#!/usr/bin/env python3
import argparse
import numpy as np

import openmc
from openmc_fusion_benchmarks import from_irdff as irdff
from openmc_fusion_benchmarks.neutron_sources import fng_source


def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batches", type=int, default=100)
    parser.add_argument("-p", "--particles", type=int, default=int(1e7))
    parser.add_argument("-s", "--threads", type=int)

    args = parser.parse_args()

    return args


def main():
    """Analysis of FNS-clean-Tungsten experiment"""

    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    ############################################################################
    # Define Materials

    # Tungsten
    mat_1 = openmc.Material(material_id=1, name='Tungsten')
    mat_1.set_density('g/cc', 18.05)
    mat_1.add_nuclide('W180', 0.0011134976037603884, 'wo')
    mat_1.add_nuclide('W182', 0.2486324290352038, 'wo')
    mat_1.add_nuclide('W183', 0.1350009119696768, 'wo')
    mat_1.add_nuclide('W184', 0.2906396609883918, 'wo')
    mat_1.add_nuclide('W186', 0.27261350040296717, 'wo')
    mat_1.add_nuclide('Ni58', 0.014111518112283878, 'wo')
    mat_1.add_nuclide('Ni60', 0.005622952208838208, 'wo')
    mat_1.add_nuclide('Ni61', 0.0002485054079123091, 'wo')
    mat_1.add_nuclide('Ni62', 0.0008053126837323057, 'wo')
    mat_1.add_nuclide('Ni64', 0.00021171158723330068, 'wo')
    mat_1.add_nuclide('Cu63', 0.021228550513624093, 'wo')
    mat_1.add_nuclide('Cu65', 0.00977144948637591, 'wo')
    mat_1.add_nuclide('Al27', 1e-08, 'wo')
    mat_1.add_nuclide('Nb93', 1e-08, 'wo')
    mat_1.add_nuclide('In113', 4.2096288100779257e-10, 'wo')
    mat_1.add_nuclide('In115', 9.579037118992209e-09, 'wo')
    mat_1.add_nuclide('Au197', 1e-08, 'wo')
    # air
    mat_2 = openmc.Material(material_id=2, name='air')
    mat_2.set_density('atom/b-cm', 4.921e-05)
    mat_2.add_nuclide('N14', 0.7886608412924202, 'ao')
    mat_2.add_nuclide('O16', 0.21133915870757974, 'ao')
    mat_2.add_nuclide('Al27', 1e-08, 'ao')
    mat_2.add_nuclide('Nb93', 1e-08, 'ao')
    mat_2.add_nuclide('In113', 4.281e-10, 'ao')
    mat_2.add_nuclide('In115', 9.571900000000001e-09, 'ao')
    mat_2.add_nuclide('Au197', 1e-08, 'ao')
    mat_2.add_nuclide('W180', 1.1999999999999999e-11, 'ao')
    mat_2.add_nuclide('W182', 2.65e-09, 'ao')
    mat_2.add_nuclide('W183', 1.4310000000000001e-09, 'ao')
    mat_2.add_nuclide('W184', 3.064e-09, 'ao')
    mat_2.add_nuclide('W186', 2.8430000000000002e-09, 'ao')

    # create materials instance
    materials = openmc.Materials([mat_1, mat_2])

    # %%

    ############################################################################
    # Define Geometry

    # surfaces
    surf_1 = openmc.ZPlane(surface_id=1, z0=-21.0, boundary_type='vacuum')
    surf_2 = openmc.ZPlane(surface_id=2, z0=-0.05)
    surf_22 = openmc.ZCylinder(surface_id=22, x0=0.0,
                               y0=0.0, r=31.45, boundary_type='vacuum')
    surf_3 = openmc.ZPlane(surface_id=3, z0=0.05)
    surf_21 = openmc.ZCylinder(surface_id=21, x0=0.0, y0=0.0, r=2.86)
    surf_4 = openmc.ZPlane(surface_id=4, z0=7.55)
    surf_5 = openmc.ZPlane(surface_id=5, z0=7.65)
    surf_6 = openmc.ZPlane(surface_id=6, z0=10.14)
    surf_7 = openmc.ZPlane(surface_id=7, z0=22.75)
    surf_8 = openmc.ZPlane(surface_id=8, z0=22.85)
    surf_9 = openmc.ZPlane(surface_id=9, z0=25.35)
    surf_10 = openmc.ZPlane(surface_id=10, z0=37.95)
    surf_11 = openmc.ZPlane(surface_id=11, z0=38.05)
    surf_12 = openmc.ZPlane(surface_id=12, z0=40.56)
    surf_13 = openmc.ZPlane(surface_id=13, z0=50.6)
    surf_14 = openmc.ZPlane(surface_id=14, z0=52.7, boundary_type='vacuum')

    # regions
    region_1 = ((+surf_1 & -surf_2 & -surf_22) |
                (+surf_2 & -surf_3 & +surf_21 & -surf_22))
    region_2 = (+surf_2 & -surf_3 & -surf_21)
    region_3 = (+surf_3 & -surf_4 & -surf_21)
    region_4 = (+surf_3 & -surf_4 & +surf_21 & -surf_22)
    region_5 = (+surf_4 & -surf_5 & -surf_21)
    region_6 = (+surf_4 & -surf_5 & +surf_21 & -surf_22)
    region_7 = (+surf_5 & -surf_6 & -surf_21)
    region_8 = (+surf_5 & -surf_6 & +surf_21 & -surf_22)
    region_9 = (+surf_6 & -surf_7 & -surf_21)
    region_10 = (+surf_6 & -surf_7 & +surf_21 & -surf_22)
    region_11 = (+surf_7 & -surf_8 & -surf_21)
    region_12 = (+surf_7 & -surf_8 & +surf_21 & -surf_22)
    region_13 = (+surf_8 & -surf_9 & -surf_21)
    region_14 = (+surf_8 & -surf_9 & +surf_21 & -surf_22)
    region_15 = (+surf_9 & -surf_10 & -surf_21)
    region_16 = (+surf_9 & -surf_10 & +surf_21 & -surf_22)
    region_17 = (+surf_10 & -surf_11 & -surf_21)
    region_18 = (+surf_10 & -surf_11 & +surf_21 & -surf_22)
    region_19 = (+surf_11 & -surf_12 & -surf_21)
    region_20 = (+surf_11 & -surf_12 & +surf_21 & -surf_22)
    region_21 = (+surf_12 & -surf_13 & -surf_21)
    region_22 = (+surf_12 & -surf_13 & +surf_21 & -surf_22)
    region_23 = (+surf_13 & -surf_14 & -surf_21)
    region_24 = (+surf_13 & -surf_14 & +surf_21 & -surf_22)
    region_25 = (-surf_1 | +surf_14 | +surf_22)

    # cells
    cell_1 = openmc.Cell(cell_id=1, region=region_1, fill=mat_2, name='Cell 1')
    cell_2 = openmc.Cell(cell_id=2, region=region_2, fill=mat_1, name='cell 2')
    cell_3 = openmc.Cell(cell_id=3, region=region_3, fill=mat_1, name='cell 3')
    cell_4 = openmc.Cell(cell_id=4, region=region_4, fill=mat_1, name='cell 4')
    cell_5 = openmc.Cell(cell_id=5, region=region_5, fill=mat_1, name='cell 5')
    cell_6 = openmc.Cell(cell_id=6, region=region_6, fill=mat_1, name='cell 6')
    cell_7 = openmc.Cell(cell_id=7, region=region_7, fill=mat_1, name='cell 7')
    cell_8 = openmc.Cell(cell_id=8, region=region_8, fill=mat_1, name='cell 8')
    cell_9 = openmc.Cell(cell_id=9, region=region_9, fill=mat_1, name='cell 9')
    cell_10 = openmc.Cell(cell_id=10, region=region_10,
                          fill=mat_1, name='cell 10')
    cell_11 = openmc.Cell(cell_id=11, region=region_11,
                          fill=mat_1, name='cell 11')
    cell_12 = openmc.Cell(cell_id=12, region=region_12,
                          fill=mat_1, name='cell 12')
    cell_13 = openmc.Cell(cell_id=13, region=region_13,
                          fill=mat_1, name='cell 13')
    cell_14 = openmc.Cell(cell_id=14, region=region_14,
                          fill=mat_1, name='cell 14')
    cell_15 = openmc.Cell(cell_id=15, region=region_15,
                          fill=mat_1, name='cell 15')
    cell_16 = openmc.Cell(cell_id=16, region=region_16,
                          fill=mat_1, name='cell 16')
    cell_17 = openmc.Cell(cell_id=17, region=region_17,
                          fill=mat_1, name='cell 17')
    cell_18 = openmc.Cell(cell_id=18, region=region_18,
                          fill=mat_1, name='cell 18')
    cell_19 = openmc.Cell(cell_id=19, region=region_19,
                          fill=mat_1, name='cell 19')
    cell_20 = openmc.Cell(cell_id=20, region=region_20,
                          fill=mat_1, name='cell 20')
    cell_21 = openmc.Cell(cell_id=21, region=region_21,
                          fill=mat_1, name='cell 21')
    cell_22 = openmc.Cell(cell_id=22, region=region_22,
                          fill=mat_1, name='cell 22')
    cell_23 = openmc.Cell(cell_id=23, region=region_23,
                          fill=mat_2, name='cell 23')
    cell_24 = openmc.Cell(cell_id=24, region=region_24,
                          fill=mat_2, name='cell 24')
    cell_25 = openmc.Cell(cell_id=25, region=region_25,
                          fill=None, name='cell 25')

    # create root universe
    universe = openmc.Universe(cells=[cell_1, cell_2, cell_3, cell_4, cell_5, cell_6, cell_7,
                                      cell_8, cell_9, cell_10, cell_11, cell_12, cell_13, cell_14,
                                      cell_15, cell_16, cell_17, cell_18, cell_19, cell_20,
                                      cell_21, cell_22, cell_23, cell_24, cell_25])

    # create geometry instance
    model.geometry = openmc.Geometry(universe)

    ############################################################################
    # Define Settings

    # source definition
    # fng source
    fng_center = (0, 0, -20)
    fng_uvw = (0., 0, 1)

    source = fng_source(center=fng_center,
                        reference_uvw=fng_uvw)

    # weight windows from wwinps
    ww = openmc.wwinp_to_wws("weight_windows.wwinp")

    # Indicate how many particles to run
    settings = openmc.Settings(run_mode='fixed source')
    settings.photon_transport = True
    settings.source = source
    settings.batches = args.batches
    settings.particles = args.particles
    settings.survival_biasing = False
    settings.weight_windows = ww
    settings.electron_treatment = 'ttb'
    settings.output = {'tallies': False}

    model.settings = settings

    ############################################################################
    # Specify Tallies

    model.tallies = openmc.Tallies()

    # energy bins for spectrum tallies
    ne213_nebins = np.array([1015913., 1068000., 1123000.,
                             1180500., 1241000., 1304500.,
                             1371500., 1442000., 1515500.,
                             1593000., 1675000., 1761000.,
                             1851000., 1946000., 2046000.,
                             2150500., 2260500., 2376500.,
                             2498500., 2626500., 2761000.,
                             2902500., 3051500., 3208000.,
                             3372500., 3545500., 3727500.,
                             3918500., 4119000., 4330500.,
                             4552500., 4786000., 5031500.,
                             5289500., 5561000., 5846000.,
                             6145500., 6460500., 6791500.,
                             7139500., 7506000., 7891000.,
                             8295500., 8721000., 9168000.,
                             9638000., 10134500., 10655000.,
                             11200000., 11775000., 12375000.,
                             13010000., 13680000., 14380000.,
                             15115000., 15890000., 16705000.,
                             17560000., 18460320.])
    prc_nebins = 1e6*np.array([3.11e-3, 3.24e-3, 3.38e-3, 3.53e-3, 3.68e-3, 3.85e-3,
                               4.02e-3, 4.20e-3, 4.38e-3, 4.58e-3, 4.78e-3, 5.00e-3,
                               5.23e-3, 5.46e-3, 5.71e-3, 5.97e-3, 6.24e-3, 6.53e-3,
                               6.83e-3, 7.14e-3, 7.47e-3, 7.81e-3, 8.17e-3, 8.55e-3,
                               8.94e-3, 9.36e-3, 9.79e-3, 1.03e-2, 1.07e-2, 1.12e-2,
                               1.18e-2, 1.23e-2, 1.29e-2, 1.35e-2, 1.41e-2, 1.48e-2,
                               1.55e-2, 1.62e-2, 1.70e-2, 1.78e-2, 1.86e-2, 1.95e-2,
                               2.04e-2, 2.14e-2, 2.24e-2, 2.34e-2, 2.45e-2, 2.57e-2,
                               2.69e-2, 2.82e-2, 2.95e-2, 3.09e-2, 3.24e-2, 3.40e-2,
                               3.56e-2, 3.73e-2, 3.90e-2, 4.09e-2, 4.29e-2, 4.49e-2,
                               4.70e-2, 4.93e-2, 5.16e-2, 5.41e-2, 5.67e-2, 5.94e-2,
                               6.22e-2, 6.52e-2, 6.83e-2, 7.16e-2, 7.50e-2, 7.86e-2,
                               8.23e-2, 8.63e-2, 9.04e-2, 9.47e-2, 9.92e-2, 1.04e-1,
                               1.09e-1, 1.14e-1, 1.20e-1, 1.25e-1, 1.31e-1, 1.38e-1,
                               1.44e-1, 1.53e-1, 1.60e-1, 1.68e-1, 1.76e-1, 1.84e-1,
                               1.93e-1, 2.02e-1, 2.12e-1, 2.22e-1, 2.33e-1, 2.44e-1,
                               2.56e-1, 2.68e-1, 2.81e-1, 2.94e-1, 3.09e-1, 3.23e-1,
                               3.39e-1, 3.55e-1, 3.72e-1, 3.90e-1, 4.09e-1, 4.28e-1,
                               4.49e-1, 4.70e-1, 4.93e-1, 5.17e-1, 5.41e-1, 5.67e-1,
                               5.95e-1, 6.23e-1, 6.53e-1, 6.84e-1, 7.17e-1, 7.51e-1,
                               7.88e-1, 8.25e-1, 8.65e-1, 9.06e-1, 9.50e-1, 9.95e-1,])
    gebins = np.array([294025., 309100., 323650., 338900.,
                       354900., 371650., 389150., 407500.,
                       426700., 446800., 467850., 489900.,
                       513000., 537150., 562450., 589000.,
                       616800., 645850., 676250., 708100.,
                       741500., 776450., 813050., 851400.,
                       891500., 933500., 977500., 1023500.,
                       1072000., 1122500., 1175000., 1230500.,
                       1288500., 1349000., 1412500., 1479500.,
                       1549500., 1622500., 1699000., 1779000.,
                       1863000., 1950500., 2042000., 2138500.,
                       2239500., 2345000., 2455500., 2571000.,
                       2692000., 2819000., 2952000., 3091000.,
                       3236500., 3389000., 3549000., 3716500.,
                       3891500., 4075000., 4267000., 4468000.,
                       4678500., 4899000., 5130000., 5371500.,
                       5624500., 5890000., 6168000., 6458500.,
                       6762500., 7081000., 7415000., 7764500.,
                       8130500., 8514000., 8915000., 9335000.,
                       9775000., 10235000., 10720000., 11225000.,
                       11750000, 12305000., 12885000., 13490000.,
                       14125000., 14795000., 15553555.])

    # filters
    # cell filters
    foil_cell_filter = openmc.CellFilter(
        [cell_2, cell_5, cell_11, cell_17, cell_23])
    spectrometer_cell_filter = [openmc.CellFilter(
        [i]) for i in [cell_5, cell_11, cell_17]]
    heatdetector_cell_filter = openmc.CellFilter(
        [cell_2, cell_5, cell_11, cell_17])
    # particle filters
    neutron_filter = openmc.ParticleFilter(['neutron'])
    photon_filter = openmc.ParticleFilter(['photon'])
    particle_filter = openmc.ParticleFilter(
        ['neutron', 'photon', 'electron', 'positron'])

    # dosimetry tallies from IRDFF-II nuclear data library
    nb93_n2n_acef = irdff.path + "dos-irdff2-4125.acef"
    al27_na_acef = irdff.path + "dos-irdff2-1325.acef"
    in115_nn_acef = irdff.path + "dos-irdff2-4931.acef"
    au197_ng_acef = irdff.path + "dos-irdff2-7925_modified.acef"
    w186_ng_acef = irdff.path + "dos-irdff2-7443_modified.acef"
    irdff_xs = [nb93_n2n_acef, al27_na_acef,
                in115_nn_acef, au197_ng_acef, w186_ng_acef]
    reactions = [11016, 107, 11004, 102, 102]
    nuclides = ['nb93', 'al27', 'in115', 'au197', 'w186']

    # create foils reaction rate cell tally
    for n, r, x in zip(nuclides, reactions, irdff_xs):
        tally1 = openmc.Tally(name=f"rr_{n}")
        irdff_xs = irdff.cross_section(x)
        multiplier = openmc.EnergyFunctionFilter.from_tabulated1d(irdff_xs[r])
        tally1.filters = [foil_cell_filter, neutron_filter, multiplier]
        tally1.scores = ["flux"]
        model.tallies.extend([tally1])

    # create neutron energy spectrum cell tally
    detector_list = ['1', '2', '3']
    for s, d in zip(spectrometer_cell_filter, detector_list):
        energy_n = openmc.EnergyFilter(ne213_nebins)
        tally2 = openmc.Tally(name=f"nspectrum_ne213_{d}")
        tally2.filters = [s, energy_n, neutron_filter]
        tally2.scores = ['flux']
        model.tallies.extend([tally2])

    for s, d in zip(spectrometer_cell_filter, detector_list):
        energy_n = openmc.EnergyFilter(prc_nebins)
        tally3 = openmc.Tally(name=f"nspectrum_prc_{d}")
        tally3.filters = [s, energy_n, neutron_filter]
        tally3.scores = ['flux']
        model.tallies.extend([tally3])

    # create gamma energy spectrum cell tally
    for s, d in zip(spectrometer_cell_filter, detector_list):
        energy_g = openmc.EnergyFilter(gebins)
        tally4 = openmc.Tally(name=f"gspectrum_bc537_{d}")
        tally4.filters = [s, energy_g, photon_filter]
        tally4.scores = ['flux']
        model.tallies.extend([tally4])

    # Tally gamma heating
    tally5 = openmc.Tally(name="nuclear_heating")
    tally5.filters = [heatdetector_cell_filter, particle_filter]
    tally5.scores = ['heating']
    model.tallies.extend([tally5])

    # define the folder names for storing the statepoints
    cwd = 'results'
    model.settings = settings

    return model.run(cwd=cwd, threads=args.threads)


if __name__ == "__main__":
    main()
