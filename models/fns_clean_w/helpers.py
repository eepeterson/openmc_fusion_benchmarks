import numpy as np

# detector volumes
detector_volume = np.ones(5) * .1*np.pi*2.86**2


ev2gy = 1.60217733e-16  # dose conversion factor from eV and Gy
# include material density for heating
_densities = np.ones(4) * 18.05

# lists
foil_list = ['nb93', 'al27', 'in115', 'au197', 'w186']
detector_list = ['1', '2', '3']


# functions
def postprocess_openmc_heating(tally_dataframe):
    """Postprocesses the openmc results for the heating simulation type

    Args:
        tally_dataframe (dataframe): dataframe with the tally results from openmc
        qtld_coeffs_table (dataframe): dataframe with qtld coefficients (has to have 'Ce', 'Cn' and 'Cp' columns)
        necessary to get results consistend with measured results for comparison

    Returns:
        dataframe: heating results postprocessed (including cell volume normalization, 
        material density normalization eV to Gy conversion and qtld coefficients implementation)
    """

    ce = 1
    cn = 1
    cp = 1

    # extract heating by particle type
    neutrons = tally_dataframe.loc[tally_dataframe['particle'] == 'neutron']
    photons = tally_dataframe.loc[tally_dataframe['particle'] == 'photon']
    electrons = tally_dataframe.loc[tally_dataframe['particle'] == 'electron']
    positrons = tally_dataframe.loc[tally_dataframe['particle'] == 'positron']
    # compute mean and std. dev.
    qn_mean = np.array(neutrons['mean'])
    qn_stddev = np.array(neutrons['std. dev.'])
    qp_mean = np.array(
        photons['mean'])+np.array(electrons['mean'])+np.array(positrons['mean'])
    qp_stddev = np.array(
        photons['std. dev.'])+np.array(electrons['std. dev.'])+np.array(positrons['std. dev.'])

    omc_qtld_mean = qn_mean*cn*ce + qp_mean*cp
    omc_qtld_mean *= ev2gy/_densities/detector_volume[:-1]
    omc_qtld_stddev = qn_stddev*cn*ce + qp_stddev*cp
    omc_qtld_stddev *= ev2gy/_densities/detector_volume[:-1]

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld_mean
    tally_dataframe['std. dev.'] = omc_qtld_stddev

    return tally_dataframe
