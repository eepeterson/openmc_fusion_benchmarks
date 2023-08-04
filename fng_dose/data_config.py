import openmc

# Set cross sections
openmc.config['cross_sections'] = '/home/peterson/xs_data/endfb80_hdf5/cross_sections.xml'

# Set chain for depletion / decay source generation
openmc.config['chain_file'] = 'chain_reduced.xml'
