from mcnptools import Mctal
import numpy as np
import h5py
import matplotlib.pyplot as plt

filename = 'mcnp_n_edit709.mcnp.m'
mctal = Mctal(filename)

tally = mctal.GetTally(4)
print(f'# of FBins: {len(tally.GetFBins())}')
print(f'# of DBins: {len(tally.GetDBins())}')
print(f'# of UBins: {len(tally.GetUBins())}')
print(f'# of SBins: {len(tally.GetSBins())}')
print(f'# of MBins: {len(tally.GetMBins())}')
print(f'# of CBins: {len(tally.GetCBins())}')
print(f'# of EBins: {len(tally.GetEBins())}')
print(f'# of TBins: {len(tally.GetTBins())}')

#print(f'FBins: {tally.GetFBins()}')
#print(f'EBins: {tally.GetEBins()}')


ebins = 1e6*np.asarray(tally.GetEBins())
lethargy = np.log(ebins[1:]/ebins[:-1])
n_ebins = len(ebins)
flux_map = {}
flux_error_map = {}
for i, cell_num in enumerate(tally.GetFBins()):
    #print(f"{i}: Cell # {int(cell_num)}")
    fluxes = np.zeros(n_ebins)
    flux_errors = np.zeros(n_ebins)
    for j in range(n_ebins):
        fluxes[j] = tally.GetValue(i, -1, -1, -1, -1, -1, j, -1)
        flux_errors[j] = tally.GetError(i, -1, -1, -1, -1, -1, j, -1)
    flux_map[int(cell_num)] = fluxes[1:]
    flux_error_map[int(cell_num)] = flux_errors[1:]

with h5py.File('mcnp_results.h5', 'w') as fh:
    mcnp_grp = fh.create_group('mcnp')
    mcnp_grp.create_dataset('ebin_edges', data=ebins)
    for cell_id in flux_map.keys():
        tmp_grp = mcnp_grp.create_group(f"cell_{cell_id}")
        tmp_grp.create_dataset('flux', data=flux_map[cell_id])
        tmp_grp.create_dataset('rel_err', data=flux_error_map[cell_id])

#cell_id = 367
cell_id = 228
#cell_id = 160
f = flux_map[cell_id] / lethargy
fsd = flux_error_map[cell_id]

fig, ax = plt.subplots()
ax.loglog(ebins[:-1], f, drawstyle='steps-post')
ax.fill_between(ebins[:-1], (1 - 1*fsd)*f, (1 + 1*fsd)*f, alpha=0.3, color='C0',
                step='post')
ax.set_xlim(1e-5, 1e8)
#ax.set_ylim(1e3, 1e14)
ax.set_xticks(np.logspace(-5,8,14))
#ax.set_yticks(np.logspace(3,14,12))
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Neutron Flux per unit lethargy (n/cm$^2$-src)')
#ax.text(1e0, 1e4, f'Total Flux: {np.sum(fluxes):1.2e}\nFast Flux (>100keV): {np.sum(fluxes[idx100:]):1.2e}', bbox=dict(boxstyle='round', fc='w'))
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3, color='k')
#ax.set_title(title)
#fig.savefig(f'{title}.png', dpi=250)
plt.show()

