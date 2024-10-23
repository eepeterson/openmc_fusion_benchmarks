import openmc


def jade_sphere(material:openmc.Material, particles:int=int(1e5), 
                photon_transport:bool=False, nreactions=None, greactions=None):

    materials = openmc.Materials([material])

    # geometry
    s1  = openmc.Sphere(r=5.)
    s2  = openmc.Sphere(r= 50.)
    s3  = openmc.Sphere(r= 60.0, boundary_type='vacuum')

    region1 = -s1
    region2 = +s1 & -s2
    region3 = +s2 & -s3

    cell1 = openmc.Cell(region=region1, fill=None)
    cell2 = openmc.Cell(region=region2, fill=material)
    cell3 = openmc.Cell(region=region3, fill=None)

    geometry = openmc.Geometry([cell1, cell2, cell3])

    # settings
    space = openmc.stats.Point()
    angle = openmc.stats.Isotropic()
    energy = openmc.stats.Discrete([14.1e6], [1])
    source = openmc.IndependentSource(space=space, angle=angle, energy=energy)

    settings = openmc.Settings(run_mode='fixed source')
    settings.particles = particles
    settings.source = source
    settings.batches = 100
    settings.output = {'tallies': False}

    tallies_file = openmc.Tallies()
    # tallies
    cell2_filter = openmc.CellFilter(cell2)
    s2_filter = openmc.SurfaceFilter(s2)
    neutron_filter = openmc.ParticleFilter(['neutron'])

    # neutron tallies
    t = openmc.Tally(name='n_flux')
    t.filters = [cell2_filter, neutron_filter]
    t.scores = ['flux']
    tallies_file.append(t)

    t = openmc.Tally(name='n_leakage')
    t.filters = [s2_filter, neutron_filter]
    t.scores = ['current']
    tallies_file.append(t)

    if nreactions is not None:
        for reaction in nreactions:
            t = openmc.Tally(name=f'n_{reaction}')
            t.filters = [cell2_filter, neutron_filter]
            t.scores = [reaction]
            tallies_file.append(t)

    if photon_transport is True:
        settings.photon_transport = True
        photon_filter = openmc.ParticleFilter(['photon'])
        t = openmc.Tally(name='g_flux')
        t.filters = [cell2_filter, photon_filter]
        t.scores = ['flux']
        tallies_file.append(t)
        t = openmc.Tally(name='g_leakage')
        t.filters = [s2_filter, photon_filter]
        t.scores = ['current']
        tallies_file.append(t)
        if nreactions is not None:
            for reaction in greactions:
                t = openmc.Tally(name=f'g_{reaction}')
                t.filters = [cell2_filter, photon_filter]
                t.scores = [reaction]
                tallies_file.append(t)
    
    return openmc.model.Model(geometry, materials, settings, tallies_file)