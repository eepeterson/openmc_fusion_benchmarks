import numpy as np

# conversions
ev2gy = 1.60217733e-16  # dose conversion factor from eV and Gy

# foil cell volumes per group
_foil_volume = .1 * 1.8**2/4 * np.pi
_tld_volume = 4/3 * np.pi * .8**3
volumes_onaxis1 = np.concatenate((np.ones(5),
                                  np.ones(4)*2,
                                  np.ones(4)*3)) * _foil_volume
volumes_onaxis2 = np.ones(11) * _foil_volume
volumes_offaxis = np.ones(15) * _foil_volume
volumes_heating = np.concatenate((np.ones(4)*_foil_volume, 
                                 np.ones(4)*2*_foil_volume, 
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

    # postprocess results
    omc_qn = []
    omc_qp = []
    for c in _cells_heating:
        c_heating = tally_dataframe.loc[tally_dataframe['cell'] == c]

        # define the right model cell volume (might be a double or a triple foil or even a sphere)
        if c in (239, 262, 285, 308):
            vol = _foil_volume * 2
        elif c in (331, 363, 386, 398):
            vol = _foil_volume * 3
        elif c in (500, 507, 514, 521):
            vol = 4/3 * np.pi * .8**3

        # cells for heating tallies are filled with either aisi316 or copper
        if c in (507, 521):
            density = 8.94  # g/cm3 - copper
        else:
            density = 7.89  # g/cm3 - aisi316

        # normalize and convert values
        # extract neutron and photons heating
        c_heating_mean = np.array(c_heating['mean']) / vol / density * ev2gy
        n_heating_mean = c_heating_mean[0]
        p_heating_mean = sum(c_heating_mean[1:])

        c_heating_stddev = np.array(
            c_heating['std. dev.']) / vol / density * ev2gy
        n_heating_stddev = c_heating_stddev[0]
        p_heating_stddev = sum(c_heating_stddev[1:])

        omc_qn.append([n_heating_mean, n_heating_stddev])
        omc_qp.append([p_heating_mean, p_heating_stddev])

    # reshape neutrons and photon heating results (both mean and std. dev.)
    omc_qn = np.array(omc_qn).T.reshape(2, 12)
    omc_qp = np.array(omc_qp).T.reshape(2, 12)

    # implement the qtld coefficients
    omc_qtld = [omc_qn[i] * cn * ce + qp * cp for i, qp in enumerate(omc_qp)]

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld[0]
    tally_dataframe['std. dev.'] = omc_qtld[1]

    return tally_dataframe
