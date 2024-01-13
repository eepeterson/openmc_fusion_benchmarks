import numpy as np

# helpers to postprocess openmc tallies

# conversions
ev2gy = 1.60217733e-16  # dose conversion factor from eV and Gy

# list of foils
foil_list = ['nb93', 'ni58_n2n', 'zr90', 'al27',
             'fe56', 'ni58_np', 'in115', 'au197', 'mn55']

# foil volumes
_fv = np.pi * .9**2
foil_volumes = {'nb93': _fv*.1, 'ni58_n2n': _fv*.1, 'zr90': _fv*.1, 'al27': _fv*.1,
                'fe56': _fv*.1, 'ni58_np': _fv*.1, 'in115': _fv*.1, 'au197': _fv*.02, 'mn55': _fv*.02}

# foil in-shield depths
_xaxis_rr1 = ['4.94', '14.94', '24.94', '34.94']
_xaxis_rr2 = ['5.04', '15.04', '25.04', '35.04']
_xaxis_rr3 = ['5.09', '15.09', '25.09', '35.09']
_xaxis_rr4 = ['5.10', '15.10', '25.10', '35.10']
xaxis_rr = {'nb93': _xaxis_rr1, 'ni58_n2n': _xaxis_rr2, 'zr90': _xaxis_rr1, 'al27': _xaxis_rr2,
            'fe56': _xaxis_rr1, 'ni58_np': _xaxis_rr2, 'in115': _xaxis_rr2, 'au197': _xaxis_rr3, 'mn55': _xaxis_rr4}

# tld - nuclear heating xaxis
xaxis_heating = ['4.96', '14.96', '24.96', '34.96']
volumes_heating = _fv * .1
_densities = 3.18  # g/cm3


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

    # # extract qtld cs, ce, cp coefficients
    # qtld_coeffs = qtld_coeffs_table
    # ce = np.array(qtld_coeffs['Ce'])
    # cn = np.array(qtld_coeffs['Cn'])
    # cp = np.array(qtld_coeffs['Cp'])

    ce = 0.1  # order of magnitude according to FNG-str TLD detector coefficient
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
    omc_qtld_mean *= ev2gy/_densities/volumes_heating
    omc_qtld_stddev = qn_stddev*cn*ce + qp_stddev*cp
    omc_qtld_stddev *= ev2gy/_densities/volumes_heating

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld_mean
    tally_dataframe['std. dev.'] = omc_qtld_stddev

    return tally_dataframe
