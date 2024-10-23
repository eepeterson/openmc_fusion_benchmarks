import openmc


def zaid_to_zam(zaid) -> tuple:
    """Converts a ZAID to Z, A, and M.

    Parameters
    ----------
    zaid : nuclide ZAID

    Returns
    -------
    tuple
        Z, A, and M values
    """
    zaid_str = str(zaid)

    # Extract Z, A, and M based on ZAID length
    if len(zaid_str) == 4:  # Handles cases like H-1 (1001)
        Z = int(zaid_str[0])  # First digit for Z
        A = int(zaid_str[1:])  # Last 3 digits for A
        M = 0  # Assuming ground state if no additional digit
    elif len(zaid_str) == 5:  # Typical ZAID with 5 digits
        Z = int(zaid_str[:2])  # First 2 digits for Z
        A = int(zaid_str[2:])  # Last 3 digits for A
        M = 0  # Assuming ground state
    elif len(zaid_str) == 6:  # ZAID with 6 digits, includes isomeric state
        Z = int(zaid_str[:3])  # First 3 digits for Z
        A = int(zaid_str[3:5])  # Next 2 digits for A
        M = int(zaid_str[5])  # Last digit represents isomeric state

    return (Z, A, M)


def get_nuclide_zaid(nuclide):
    """Gets the ZAID of a nuclide from its name as a GNDS string (e.g. 'H1',
    'U238) or as a ZAM tuple (e.g. (1, 1, 0)).
    See openmc.data.zam() for more details.

    Parameters
    ----------
    nuclide : int or str or tuple
        The nuclide identifier or name

    Returns
    -------
    int
        The ZAID of the nuclide
    """
    if type(nuclide) == int:
        return nuclide
    elif type(nuclide) == str:
        return openmc.data.zam(nuclide)[0]*1000 + openmc.data.zam(nuclide)[1]
    elif type(nuclide) == tuple:
        return nuclide[0]*1000 + nuclide[1]


def get_nuclide_gnds(nuclide):
    """Gets the GNDS name from a nuclide ZAID

    Parameters
    ----------
    nuclide : str or int
        The nuclide ZAID or GNDS name

    Returns
    -------
    str
        The GNDS name of the nuclide
    """
    if type(nuclide) == str:
        return nuclide
    elif type(nuclide) == int:
        zam = zaid_to_zam(nuclide)
        return openmc.data.gnds_name(zam[0], zam[1], zam[2])
