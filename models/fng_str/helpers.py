import numpy as np

# conversions
ev2gy = 1.60217733e-16  # dose conversion factor from eV and Gy

# foil cell volumes per group
_foil_volume = .1 * .9**2 * np.pi
_tld_volume = 4/3 * np.pi * .8**3
volumes_onaxis1 = np.concatenate((np.ones(5),
                                  np.ones(4)*2,
                                  np.ones(4)*3)) * _foil_volume
volumes_onaxis2 = np.ones(11) * _foil_volume
volumes_offaxis = np.ones(15) * _foil_volume
volumes_heating = np.concatenate((np.ones(4)*2*_foil_volume, 
                                 np.ones(4)*3*_foil_volume, 
                                 np.ones(4)*_tld_volume))
# include material density for heating
_densities = np.ones(12) * 7.89
_densities[-1] = _densities[-3] = 8.94

# heating simulation data
_cells_heating = [239, 262, 285, 308, 331, 363, 386, 398, 500, 507, 514, 521]


def postprocess_openmc_heating(tally_dataframe, qtld_coeffs_table):
    """Postprocesses the openmc results for the heating simulation type

    Args:
        tally_dataframe (dataframe): dataframe with the tally results from openmc
        qtld_coeffs_table (dataframe): dataframe with qtld coefficients (has to have 'Ce', 'Cn' and 'Cp' columns)
        necessary to get results consistend with measured results for comparison

    Returns:
        dataframe: heating results postprocessed (including cell volume normalization, 
        material density normalization eV to Gy conversion and qtld coefficients implementation)
    """

    # extract qtld cs, ce, cp coefficients
    qtld_coeffs = qtld_coeffs_table
    ce = np.array(qtld_coeffs['Ce'])
    cn = np.array(qtld_coeffs['Cn'])
    cp = np.array(qtld_coeffs['Cp'])

    # extract heating by particle type
    neutrons = tally_dataframe.loc[tally_dataframe['particle'] == 'neutron']
    photons = tally_dataframe.loc[tally_dataframe['particle'] == 'photon']
    electrons = tally_dataframe.loc[tally_dataframe['particle'] == 'electron']
    positrons = tally_dataframe.loc[tally_dataframe['particle'] == 'positron']
    # compute mean and std. dev.
    qn_mean = np.array(neutrons['mean'])
    qn_stddev = np.array(neutrons['std. dev.'])
    qp_mean = np.array(photons['mean'])+np.array(electrons['mean'])+np.array(positrons['mean'])
    qp_stddev = np.array(photons['std. dev.'])+np.array(electrons['std. dev.'])+np.array(positrons['std. dev.'])

    omc_qtld_mean = qn_mean*cn*ce + qp_mean*cp
    omc_qtld_mean *= ev2gy/_densities/volumes_heating
    omc_qtld_stddev = qn_stddev*cn*ce + qp_stddev*cp
    omc_qtld_stddev *= ev2gy/_densities/volumes_heating

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld_mean
    tally_dataframe['std. dev.'] = omc_qtld_stddev

    return tally_dataframe
