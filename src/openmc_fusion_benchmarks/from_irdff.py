import openmc.data

path = "../../src/openmc_fusion_benchmarks/data/irdff2_xs/"


def cross_section(irdff_file_path: str):
    """Generates cross section data from IRDFF-II files
    from this discussion ad related notebook:
    https://openmc.discourse.group/t/using-irdff-ii-cross-section-data-in-openmc/1950

    Parameters
    ----------
    irdff_file_path : str
        takes in a string with the path and the name of the .acef file
        specific for the nuclide in IRDFF-II nuclear data library

    Returns
    -------
    openmc.data.Tabulated1D
        IRDFF-II tabulated cross section data for a given nuclide
        and reaction
    """
    ace_table = openmc.data.ace.get_table(irdff_file_path)
    nxs = ace_table.nxs
    jxs = ace_table.jxs
    xss = ace_table.xss

    # Get MT values and locators for cross section blocks
    lmt = jxs[3]
    nmt = nxs[4]
    lxs = jxs[6]
    mts = xss[lmt: lmt + nmt].astype(int)
    locators = xss[lxs: lxs + nmt].astype(int)

    # Create dictionary mapping MT to Tabulated1D object
    cross_sections = {}
    for mt, loca in zip(mts, locators):
        # Determine starting index on energy grid
        nr = int(xss[jxs[7] + loca - 1])
        if nr == 0:
            breakpoints = None
            interpolation = None
        else:
            breakpoints = xss[jxs[7] + loca: jxs[7] + loca + nr].astype(int)
            interpolation = xss[jxs[7] + loca +
                                nr: jxs[7] + loca + 2 * nr].astype(int)

        # Determine number of energies in reaction
        ne = int(xss[jxs[7] + loca + 2 * nr])

        # Read reaction cross section
        start = jxs[7] + loca + 1 + 2 * nr
        energy = xss[start: start + ne] * 1e6
        xs = xss[start + ne: start + 2 * ne]

        cross_sections[mt] = openmc.data.Tabulated1D(
            energy, xs, breakpoints, interpolation
        )

    return cross_sections
