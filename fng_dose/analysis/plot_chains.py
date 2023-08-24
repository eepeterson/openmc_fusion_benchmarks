import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle, Patch
import openmc.deplete
import openmc.data

# customizations
mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
mpl.rcParams['axes.axisbelow'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['font.size'] = 12.0

chain_endf80_reduced = openmc.deplete.Chain.from_xml('../chain_endfb80_reduced.xml')
chain_endf71_reduced = openmc.deplete.Chain.from_xml('../chain_reduced.xml')

def get_patch(z, a, **kwargs):
    n = a - z
    return Rectangle((n-0.5, z-0.5), 1, 1, linewidth=0.001, **kwargs)

fig, ax = plt.subplots()

style_endf71 = {'color': '#bcd6e9'}
style_endf80 = {'color': '#1f77b4'}
style_casl = {'color': 'orange'}

for nuc in chain_endf80_reduced.nuclides:
    z, a, m = openmc.data.zam(nuc.name)
    if m > 0: continue
    if nuc.name not in chain_endf71_reduced:
        ax.add_patch(get_patch(z, a, **style_endf80))

for nuc in chain_endf71_reduced.nuclides:
    z, a, m = openmc.data.zam(nuc.name)
    if m > 0: continue
    ax.add_patch(get_patch(z, a, **style_endf71))

ax.set_xlim(0, 140)
ax.set_ylim(0, 100)
ax.set_xlabel('N')
ax.set_ylabel('Z')
ax.grid()
legend_elements = [
    Patch(label='ENDF/B-VII.1', **style_endf71),
    Patch(label='ENDF/B-VIII.0', **style_endf80),
#    Patch(label='CASL', **style_casl)
]
ax.legend(handles=legend_elements, framealpha=1.0)

plt.savefig('activation_chains.pdf')
