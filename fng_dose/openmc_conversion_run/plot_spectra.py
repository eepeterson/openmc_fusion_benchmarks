import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams["mathtext.default"] = "regular"

# Maps of Cell ID to tuple of (flux, rel_err)
mcnp_fluxes = {}
openmc_fluxes = {}
rms_map = {}

with h5py.File('flux_spectrum_results.h5', 'r') as fh:
    ebins = fh['openmc/ebin_edges'][:]
    lethargy = np.log(ebins[1:] / ebins[:-1])
    for cell_tag in fh['openmc'].keys():
        if "cell" not in cell_tag:
            continue
        cell_id = int(cell_tag[5:])
        mcnp_f = fh[f"mcnp/{cell_tag}/flux"][:]
        mcnp_re = fh[f"mcnp/{cell_tag}/rel_err"][:]
        mcnp_fluxes[cell_id] = (mcnp_f, mcnp_re)
        openmc_f = fh[f"openmc/{cell_tag}/flux"][:]
        openmc_re = fh[f"openmc/{cell_tag}/rel_err"][:]
        openmc_re[~np.isfinite(openmc_re)] = 0.0
        openmc_fluxes[cell_id] = (openmc_f, openmc_re)

for cell_id in openmc_fluxes.keys():
    f1 = openmc_fluxes[cell_id][0]
    f2 = mcnp_fluxes[cell_id][0]
    rms_map[cell_id] = (f1 - f2)**2 / np.sum(f2)**2



fig, ax = plt.subplots()
for cell_id in openmc_fluxes.keys():
    f1 = openmc_fluxes[cell_id][0]
    fsd1 = openmc_fluxes[cell_id][1]
    f2 = mcnp_fluxes[cell_id][0]
    fsd2 = mcnp_fluxes[cell_id][1]
    nz_idx = np.where((fsd1 > 0) & (fsd1 < 0.5) & (fsd2 > 0) & (fsd2 < 0.5))
    c_over_e = np.abs((f1 - f2)) / f2
    sorted_ce = sorted(c_over_e[nz_idx])
    pct_ce = np.linspace(0, 1, len(sorted_ce))
    ax.loglog(sorted_ce, pct_ce*100, color='C0', alpha=0.3)

ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.3, color='k')
ax.grid(True, linestyle='--', which='minor', linewidth=0.5, alpha=0.15, color='k')

ax.set_xlabel('Relative Difference')
ax.set_ylabel('Energy Bin Percentage (%)')
ax.set_ylim(.2, 200)
ax.set_xlim(1e-4, 2)
plt.show()

hist_vals = [np.sqrt(np.sum(val)) for val in rms_map.values()]
fig, ax = plt.subplots()
ax.hist(hist_vals, bins=30)
plt.show()

plot_cells = [160, 228, 367]
#plot_cells = openmc_fluxes.keys()

for cell_id in plot_cells:
    f1 = openmc_fluxes[cell_id][0]
    fsd1 = openmc_fluxes[cell_id][1]
    f2 = mcnp_fluxes[cell_id][0]
    fsd2 = mcnp_fluxes[cell_id][1]
    c_over_e = (f1 - f2) / f2

    fig, (ax0, ax1) = plt.subplots(2, 1, gridspec_kw={'height_ratios':(2, 1)},
                                   figsize=(9,6))
    ax0.loglog(ebins[:-1], f2, drawstyle='steps-post', label="MCNP 6.2")
    ax0.fill_between(ebins[:-1], (1 - 1*fsd2)*f2, (1 + 1*fsd2)*f2, alpha=0.3, color='C0',
                    step='post')
    ax0.loglog(ebins[:-1], f1, drawstyle='steps-post', label="OpenMC 0.13.4-dev")
    ax0.set_xticks(np.logspace(-4,8,13))
    ax0.set_xlim(1e-4, 2e7)
    ax0.get_xaxis().set_ticklabels([])
    ax0.set_ylabel('Neutron flux per unit lethargy (n/cm$^2$-src)')
    ax0.legend()
    #ax0.text(1e0, 1e4, f'Total Flux: {np.sum(fluxes):1.2e}\nFast Flux (>100keV): {np.sum(fluxes[idx100:]):1.2e}', bbox=dict(boxstyle='round', fc='w'))
    ax0.grid(True, linestyle='--', linewidth=0.5, alpha=0.3, color='k')
    title = f"Cell {cell_id}"
    ax0.set_title(title)

    ax1.fill_between(ebins[:-1], -3*fsd2, 3*fsd2, alpha=0.3, color='C0', step='post')
    ax1.fill_between(ebins[:-1], -2*fsd2, 2*fsd2, alpha=0.3, color='C0', step='post')
    ax1.fill_between(ebins[:-1], -1*fsd2, 1*fsd2, alpha=0.3, color='C0', step='post')
    ax1.semilogx(ebins[:-1], c_over_e, drawstyle='steps-post', color='C1')
    ax1.set_ylim(-.1, .1)
    ax1.set_xticks(np.logspace(-4,8,13))
    ax1.set_xlim(1e-4, 2e7)
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel(r'$\frac{OpenMC}{MCNP} - 1$')
    #ax1.text(1e0, 1e4, f'Total Flux: {np.sum(fluxes):1.2e}\nFast Flux (>100keV): {np.sum(fluxes[idx100:]):1.2e}', bbox=dict(boxstyle='round', fc='w'))
    #ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.3, color='k')
    fig.tight_layout()
    plt.show()
    #print(title)
    #fig.savefig(f'neutron_spectrum_figures/{title}.png', dpi=150)
    #plt.close()
