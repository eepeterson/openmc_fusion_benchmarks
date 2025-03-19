#!/usr/bin/env python3
import numpy as np
import argparse

import openmc
from openmc_fusion_benchmarks import from_irdff as irdff
import openmc.lib
import os
from pathlib import Path
import glob
#import pandas as pd
import openmc.mgxs as mgxs

def _parse_args():
    """Parse and return commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--batches", type=int, default=20)
    parser.add_argument("-p", "--particles", type=int, default=int(1e7))
    parser.add_argument("-s", "--threads", type=int,default=5)
    parser.add_argument("-c", "--cwd", type=str)
    parser.add_argument("-I", "--Problem_I", action='store_true',default=False)
    parser.add_argument("-w", "--weight_windows",action='store_true',default=False)

    args = parser.parse_args()

    return args


def gen_ww(tally_id):

    openmc.lib.init([])
    openmc.lib.simulation_init()

    tally = openmc.lib.tallies[tally_id]

    wws = openmc.lib.WeightWindows.from_tally(tally)

    iterations = 10
    openmc.lib.settings.particles = int(1e6)
    n_batches = openmc.lib.settings.get_batches()

    for i in range(iterations):

        openmc.lib.run() #run the simulation

        os.rename(f'statepoint.{n_batches}.h5', f'statepoint.{n_batches}.{i}.h5')
        print("Updating")
        wws.update_magic(tally, value='rel_err', threshold=1.0) #updates the weight window with the latest tally results
        print("Done updating")

        openmc.lib.settings.weight_windows_on = True #turns on weight windows to ensure they are used
        if i % 2 == 0: #this doubles the amount of particles were iteration
            openmc.lib.settings.particles *= 2
        openmc.lib.reset()

        print("Writing weight windows")
        #write out the weight window maps for plotting later
        openmc.lib.export_weight_windows(filename = f'weight_windows{i}.h5')

    openmc.lib.finalize()

def main():
    """Analysis of Kobayashi Dogleg benchmark"""

    #This function builds the geeometry & materials, and can
    # be used to generate the weight windows files needed for
    # Problems 3Ci and 3Cii in the Kobayashi Dogleg benchmark.
    # The first run-through should generate the weight windows,
    # and the second should generate the 'production' runs.

    # The production runs follow the examples of this library
    # in returning the model.run call using subprocess. However,
    # the weight windows generation was written using openmc.lib
    # so that happens here within main() and returns nothing


    # Parse commandline arguments
    args = _parse_args()

    # Instantiate Model object
    model = openmc.Model()

    # define materials

    #Make a 1-group object
    # the energy group boundaries in eV
    groups = mgxs.EnergyGroups(group_edges = np.array([0., 20.0e6]))


    ############################################################################
    # Build Universe

    # for my universe cell
    plane1 = openmc.XPlane(x0= -0.0, boundary_type= 'reflective')
    plane2 = openmc.XPlane(x0= +60, boundary_type= 'vacuum')
    plane3 = openmc.YPlane(y0= -0.0, boundary_type= 'reflective')
    plane4 = openmc.YPlane(y0= +100, boundary_type= 'vacuum')
    plane5 = openmc.ZPlane(z0= -0.0, boundary_type= 'reflective')
    plane6 = openmc.ZPlane(z0= +60, boundary_type= 'vacuum')


    box = +plane1 & -plane2 & +plane3 & -plane4 & +plane5 & -plane6

    Ucell = openmc.Cell(name='base universe', region = box)

    # domain: the domain for spatial homogenization aka
    # averaging the cross sections over volume and flux spectrum

    total = mgxs.TotalXS(domain=Ucell, energy_groups=groups)
    absorption = mgxs.AbsorptionXS(domain=Ucell, energy_groups=groups)
    scattering = mgxs.ScatterXS(domain=Ucell, energy_groups=groups)

    ebins = [1e-5, 20.0e6]
    groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

    ############################################################################
    # Build Multi-group cross-sections

    #Problem_I = no scattering (only absorption)
    #Problem_II = 50% scattering in absorber


    void_sigma_a = 0.5e-4
    void_sigma_s = 0.5e-4
    if args.Problem_I:
        void_sigma_s = 0
        void_sigma_a = 1.0e-4

    void_mat_data = openmc.XSdata('void', groups)
    void_mat_data.order = 0
    void_mat_data.set_total([void_sigma_a + void_sigma_s])
    void_mat_data.set_absorption([void_sigma_a])
    void_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[void_sigma_s]]]), 0, 3))

    absorber_sigma_a = 0.05
    absorber_sigma_s = 0.05
    if args.Problem_I:
        absorber_sigma_a = 0.1
        absorber_sigma_s = 0.0

    absorber_mat_data = openmc.XSdata('absorber', groups)
    absorber_mat_data.order = 0
    absorber_mat_data.set_total([absorber_sigma_a + absorber_sigma_s])
    absorber_mat_data.set_absorption([absorber_sigma_a])
    absorber_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[absorber_sigma_s]]]), 0, 3))

    source_sigma_a = 0.05
    source_sigma_s = 0.05
    if args.Problem_I:
        source_sigma_a = 0.1
        source_sigma_s = 0.0

    source_mat_data = openmc.XSdata('source', groups)
    source_mat_data.order = 0
    source_mat_data.set_total([source_sigma_a + source_sigma_s])
    source_mat_data.set_absorption([source_sigma_a])
    source_mat_data.set_scatter_matrix(
        np.rollaxis(np.array([[[source_sigma_s]]]), 0, 3))

    #Output the cross-sections we just created to a hdf5 file, which we
    # will later tell OpenMC to load
    mg_cross_sections_file = openmc.MGXSLibrary(groups)
    mg_cross_sections_file.add_xsdatas(
        [source_mat_data, void_mat_data, absorber_mat_data])
    if args.Problem_I:
        mg_cross_sections_file.export_to_hdf5("mgxs_Problem_I.h5")
    else:
        mg_cross_sections_file.export_to_hdf5("mgxs_Problem_II.h5")
    # Define Macroscopic Materials

    source_data = openmc.Macroscopic('source')
    void_data = openmc.Macroscopic('void')
    absorber_data = openmc.Macroscopic('absorber')

    Source = openmc.Material(name='source_mat')
    Source.set_density('macro', 1.0)
    Source.add_macroscopic('source')

    Absorber = openmc.Material(name='absorber_mat')
    Absorber.set_density('macro', 1.0)
    Absorber.add_macroscopic('absorber')

    Void = openmc.Material(name='void_mat')
    Void.set_density('macro', 1.0)
    Void.add_macroscopic('void')

    # Instantiate a Materials collection and export to XML
    materials = openmc.Materials([Source,Absorber,Void])
    if args.Problem_I:
        materials.cross_sections = "mgxs_Problem_I.h5"
    else:
        materials.cross_sections = "mgxs_Problem_II.h5"
    materials.export_to_xml()
    model.materials = materials


    ###############################################################################
    # Build Geometry

     #Produce the void dog leg regions:

    #Void Segment 1:

    plane7 = openmc.XPlane(x0= +10, boundary_type= 'transmission')
    plane12 = openmc.YPlane(y0= +10, boundary_type= 'transmission')
    plane8 = openmc.YPlane(y0= +60, boundary_type= 'transmission')
    plane9 = openmc.ZPlane(z0= +10, boundary_type= 'transmission')

    segment1 = +plane1 & -plane7 & +plane12 & -plane8 & +plane5 & -plane9 #(0,10,10,60,0,10)

    #Void Segment 2:

    plane10 = openmc.XPlane(x0= +40, boundary_type= 'transmission')
    plane11 = openmc.YPlane(y0= +50, boundary_type= 'transmission')

    segment2 = +plane7 & -plane10 & +plane11 & -plane8 & +plane5 & -plane9 #(10,40,50,60,0,10)


    #Void Segment 3:
    # can keep parallelpiped since all transmission boundaries
    segment3 = openmc.model.RectangularParallelepiped(30,40,50,60,10,40)

    #Void Segment 4:
    # can keep parallelpiped since all transmission boundaries
    segment4 = openmc.model.RectangularParallelepiped(30,40,60,99.999,30,40)

    box2 = +plane1 & -plane7 & +plane3 & -plane12 & +plane5 & -plane9 #(0,10,0,10,0,10)

    box0 = +plane1 & -plane2 & +plane3 & -plane4 & +plane5 & -plane6


    source_cell = openmc.Cell(fill = Source,
                          region = box2,
                          name = 'source region')

    void_cell = openmc.Cell(fill = Void,
                            region = (segment1 | segment2 | -segment3 | -segment4),
                            name = 'void region')

    absorber_cell = openmc.Cell(fill = Absorber,
                                region = box0 & ~void_cell.region & ~source_cell.region,
                                name = 'absorber region')


    # Fill the Cell with the Material
    Ucell.fill = Void

    # Create root universe
    root_universe = openmc.Universe(name='root universe',
                                    cells=[source_cell, void_cell, absorber_cell])

    # Create Geometry and set root Universe
    geometry = openmc.Geometry(root_universe)
    geometry.merge_surfaces=True

    ## Export to "geometry.xml"
    geometry.export_to_xml()
    model.geometry = geometry

    ###############################################################################
    # Define problem settings

    # Indicate how many particles to run
    settings = openmc.Settings(run_mode='fixed source')
    settings.energy_mode = "multi-group"
    settings.batches = args.batches
    settings.particles = args.particles
    # create an initial uniform spatial source distribution
    bounds = [0,0,0,10,10,10] # match the geometry bounds for my source [xmin,ymin,zmin,xmax,ymax,zmax]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])
    energy_mono = openmc.stats.Discrete([1.0], [1.0]) # specifying a monenergentic source
    settings.source = openmc.IndependentSource(energy = energy_mono, space = uniform_dist)
    settings.output = {'tallies': False}

    ############################################################################
    # define tallies

    model.tallies = openmc.Tallies()

    if args.weight_windows:

        geo_mesh = openmc.RegularMesh()
        geo_mesh.dimension = [60,100,60]
        geo_mesh.lower_left = [0,0,0]
        geo_mesh.upper_right =[60,100,60]

        geo_mesh_filter = openmc.MeshFilter(geo_mesh)
        flux_tally = openmc.Tally(name= 'Flux Tally')
        flux_tally.filters = [geo_mesh_filter]
        flux_tally.scores = ["flux"]

        # Set ID number to access during openmc run:
        flux_tally.id = 4
        model.tallies.append(flux_tally)
        #must export to xml so it gets picked up in the function above
        model.tallies.export_to_xml()


    else:
        a_mesh = openmc.RegularMesh()
        a_mesh.dimension = [1,10,1] # number of bins
        a_mesh.lower_left = [4.5, 4.5, 4.5] # physical "corners" of mesh
        a_mesh.upper_right =[5.5, 95.5, 5.5]

        b_mesh = openmc.RegularMesh()
        b_mesh.dimension = [6,1,1]
        b_mesh.lower_left = [4.5 , 54.5 , 4.5]
        b_mesh.upper_right =[55.5, 55.5, 5.5]

        c_mesh = openmc.RegularMesh()
        c_mesh.dimension = [6,1,1]
        c_mesh.lower_left = [4.5, 94.5, 34.5]
        c_mesh.upper_right =[55.5, 95.5 , 35.5]

        a_mesh_filter = openmc.MeshFilter(a_mesh)
        b_mesh_filter = openmc.MeshFilter(b_mesh)
        c_mesh_filter = openmc.MeshFilter(c_mesh)

        tallyname = 'Problem_II_3A'
        if args.Problem_I:
            tallyname =  'Problem_I_3A'
        a_tally = openmc.Tally(name = tallyname)
        a_tally.filters = [a_mesh_filter]
        a_tally.scores = ['flux']

        tallyname = 'Problem_II_3B'
        if args.Problem_I:
            tallyname = 'Problem_I_3B'
        b_tally = openmc.Tally(name = tallyname)
        b_tally.filters = [b_mesh_filter]
        b_tally.scores = ['flux']

        tallyname = 'Problem_II_3C'
        if  args.Problem_I :
            tallyname = 'Problem_I_3C'
        c_tally = openmc.Tally(name = tallyname)
        c_tally.filters = [c_mesh_filter]
        c_tally.scores = ['flux']

        model.tallies.extend([a_tally,b_tally,c_tally])
        #model.tallies.export_to_xml()

    model.settings = settings
    settings.export_to_xml()

    # define the folder names for storing the statepoints
    if args.Problem_I:
        cwd = 'Problem_I/'
    else:
        cwd = 'Problem_II/'

    Path(cwd).mkdir(exist_ok=True)

    if args.weight_windows:

        ############################################################################
        # Switch for generating weight windows with the Magic Method
        # or for running with the generated weight windows

        gen_ww(4) #use the tally id defined above

        for j in glob.glob("*.xml"):
            os.rename(j,cwd+j)
        for j in glob.glob("*.h5"):
            os.rename(j,cwd+j)

        return 0

    else:
        #Add the weight windows information to settings
        settings.weight_windows_file='weight_windows9.h5'
        model.settings = settings

        return model.run(cwd=cwd, threads=args.threads)



if __name__ == "__main__":
    main()
