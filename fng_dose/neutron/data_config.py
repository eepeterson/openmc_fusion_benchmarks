import openmc

# Set cross sections
openmc.config['cross_sections'] = '/opt/data/hdf5/endfb-viii.0-hdf5/cross_sections.xml'

# Set chain for depletion / decay source generation
openmc.config['chain_file'] = 'chain_reduced.xml'
