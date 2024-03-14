import openmc
import numpy as np
import os
import pandas as pd

# data necessary for postprocessing
_sim_types = ['reaction_rates_onaxis', 'reaction_rates_offaxis', 'heating']
# reaction rate simulations
# possible group types for the reaction rate simulation types
_groups = ['onaxis_group1', 'onaxis_group2', 'offaxis']
# list of foil types
_foil_list = ['nb93', 'al27', 'ni58', 'au197']
# openmc model foil cell id per group
_cells_group1 = [135, 158, 181, 204, 602,
                 239, 262, 285, 308, 331, 363, 386, 398]
_cells_group2 = [605, 606, 607, 608, 609, 610, 611, 612, 602, 603, 604]
_cells_offaxis = [135, 158, 181, 204, 605, 606,
                  607, 608, 609, 610, 611, 612, 602, 603, 604]
# foil cell volumes per group
_foil_volume = .1 * 1.8**2/4 * np.pi
_volumes_group1 = np.concatenate(
    (np.ones(5), np.ones(4)*2, np.ones(4)*3)) * _foil_volume
_volumes_group2 = np.ones(11) * _foil_volume
_volumes_offaxis = np.ones(15) * _foil_volume

# heating simulation data
_cells_heating = [239, 262, 285, 308, 331, 363, 386, 398, 500, 507, 514, 521]
_ev2gy = 1.60217733e-16  # dose conversion factor from eV and Gy
# extract qtld cs, ce, cp coefficients
_ce = np.array([0.1, 0.07, 0.05, 0.05, 0.05, 0.05,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
_cn = np.array([2.7, 2.78, 2.9, 3., 3.09, 3.19,
                3.29, 3.39, 3.5,  4.54, 3.65, 5.12])
_cp = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.05, 1., 1.05])

statepoint_file = r"statepoint.100.h5"


def get_openmc_tally(simulation_type: str, statepoint_file: str, tally_name: str):

    # raise error if wrong simulation_type
    if simulation_type not in _sim_types:
        raise ValueError(
            "Invalid simulation_type. Expected one of: %s" % _sim_types)

    # enter the right folder according to the simulation type, read the openmc statepoint file and convert to dataframe
    try:
        os.chdir(simulation_type)
        sp = openmc.StatePoint(statepoint_file)
        tally_dataframe = sp.get_tally(name=tally_name).get_pandas_dataframe()
        os.chdir('..')
        return tally_dataframe
    except FileNotFoundError:
        pass


def postprocess_openmc_foils(tally_dataframe, foil_group: str):

    if foil_group not in _groups:
        raise ValueError("Invalid foil_group. Expected one of: %s" % _groups)

    df = tally_dataframe.loc[(tally_dataframe['particle'] == 'neutron')]

    if foil_group == _groups[0]:
        df = df[df['cell'].isin(_cells_group1)]
        df['mean'] = df['mean'] / _volumes_group1
        df['std. dev.'] = df['std. dev.'] / _volumes_group1
    elif foil_group == _groups[1]:
        df = df[df['cell'].isin(_cells_group2)]
        df['mean'] = df['mean'] / _volumes_group2
        df['std. dev.'] = df['std. dev.'] / _volumes_group2
    elif foil_group == _groups[2]:
        df = df[df['cell'].isin(_cells_offaxis)]
        df['mean'] = df['mean'] / _volumes_offaxis
        df['std. dev.'] = df['std. dev.'] / _volumes_offaxis

    tally_dataframe = df[['cell', 'mean', 'std. dev.']]

    return tally_dataframe


def postprocess_openmc_heating(tally_dataframe):

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
        c_heating_mean = np.array(c_heating['mean']) / vol / density * _ev2gy
        n_heating_mean = c_heating_mean[0]
        p_heating_mean = sum(c_heating_mean[1:])

        c_heating_stddev = np.array(
            c_heating['std. dev.']) / vol / density * _ev2gy
        n_heating_stddev = c_heating_stddev[0]
        p_heating_stddev = sum(c_heating_stddev[1:])

        omc_qn.append([n_heating_mean, n_heating_stddev])
        omc_qp.append([p_heating_mean, p_heating_stddev])

    # reshape neutrons and photon heating results (both mean and std. dev.)
    omc_qn = np.array(omc_qn).T.reshape(2, 12)
    omc_qp = np.array(omc_qp).T.reshape(2, 12)

    # implement the qtld coefficients
    omc_qtld = [omc_qn[i] * _cn * _ce + qp *
                _cp for i, qp in enumerate(omc_qp)]

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld[0]
    tally_dataframe['std. dev.'] = omc_qtld[1]

    return tally_dataframe[['cell', 'mean', 'std. dev.']]


# postprocess and write csv files
try:
    for foil in _foil_list:
        tdf = get_openmc_tally('reaction_rates_onaxis',
                               'statepoint.100.h5', f'{foil}_irdff_rr')
        group1 = postprocess_openmc_foils(
            tdf, 'onaxis_group1').to_csv(f'onaxis_group1_{foil}.csv')
        group2 = postprocess_openmc_foils(
            tdf, 'onaxis_group2').to_csv(f'onaxis_group2_{foil}.csv')
        tdf = get_openmc_tally('reaction_rates_offaxis',
                               'statepoint.100.h5', f'{foil}_irdff_rr')
        offaxis = postprocess_openmc_foils(
            tdf, 'offaxis').to_csv(f'offaxis_{foil}.csv')
    tdf = get_openmc_tally('heating',
                           'statepoint.100.h5', 'heating_dose')
    heating = postprocess_openmc_heating(tdf).to_csv('heating_dose.csv')
except:
    pass
