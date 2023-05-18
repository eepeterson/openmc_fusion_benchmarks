import numpy as np
import math
import matplotlib.pyplot as plt
import os

import openmc

# data necessary for postprocessing
_sim_types = ['reaction_rates_onaxis', 'reaction_rates_offaxis', 'heating']
# reaction rate simulations
# possible group types for the reaction rate simulation types
_groups = ['onaxis_group1', 'onaxis_group2', 'offaxis']
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


# functions for postprocessing and visualization
def _get_floor_ceiling(values):
    """Finds the global minimum and maximum values that can be used to frame a plot

    Args:
        values (array): values that will be plotted

    Returns:
        float, float: global minimum and global maximum values present in the input array
    """

    min_value = min([np.nanmin(i) for i in values])
    max_value = max([np.nanmax(i) for i in values])
    min_oom = math.floor(math.log(min_value, 10))
    max_oom = math.floor(math.log(max_value, 10))
    return 10**(min_oom),  10**(max_oom+1)


def mean(dataset):
    return dataset['mean']


def std_dev(dataset):
    return dataset['std. dev.']


def rel_std_dev(dataset):
    return dataset['std. dev.'] / dataset['mean']


def mean_ratio(dataset_1, dataset_2):
    return dataset_1 / dataset_2


def _3sigma_area(ax, ticks, std_dev):
    """Generates errorbars as grey areas on existing plot

    Args:
        ax (matplotlib.pyplot.axes): plot already existing
        ticks (array): array of values desctibing the plot xaxis
        std_dev (array): values corresponding to 1 standard deviation
    """
    # plotting the three grey areas
    ax.fill_between(ticks, 1 - std_dev, 1 + std_dev, color='k', alpha=.1)
    ax.fill_between(ticks, 1 - 2*std_dev, 1 + 2*std_dev, color='k', alpha=.1)
    ax.fill_between(ticks, 1 - 3*std_dev, 1 + 3*std_dev, color='k', alpha=.1)


def _print_3sigma(ax, mean, std_dev, ticks):
    """Generates the 1,2,3simga strings near the grey areas for an existing plot

    Args:
        ax (matplotlib.pyplot.axes): plot already existing
        mean (array): values corresponding to the mean
        std_dev (array): values corresponding to 1 standard deviation
        ticks (array): array of values desctibing the plot xaxis
    """
    # get index of last valid value for measured data for placing 1,2,3sigma strings near the grey areas
    sigma_index = (~np.isnan(mean)).cumsum(0).argmax(0)
    last_sigma = std_dev[sigma_index]

    # trying to avoid interference between the 3 strings
    if last_sigma < .08:
        last_sigma = .08
        sigma_start = .93

    # string annootation in the plot
    ax.annotate('1\u03C3', [ticks[sigma_index] + 0.07,
                1. + 1*last_sigma], fontsize=10, clip_on=False)
    ax.annotate('2\u03C3', [ticks[sigma_index] + 0.07,
                1. + 2*last_sigma], fontsize=10, clip_on=False)
    ax.annotate('3\u03C3', [ticks[sigma_index] + 0.07,
                1. + 3*last_sigma], fontsize=10, clip_on=False)


class VisualizeResults:
    """Embeds all necessary functions for plotting results and C/E results for comparison
    __init__() already generates the plot istance
    """

    def __init__(self):

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6, 5),
                                       gridspec_kw={'height_ratios': [2, 1.25]}, constrained_layout=True)

        self.fig = fig
        self.ax1 = ax1
        self.ax2 = ax2

    def add_measured_data(self, measured_data, ylabel: str, dtype_label=''):
        """Includes measured data in the already existing plot generated in the class initialization

        Args:
            measured_data (dataframe or dataset): table with the measured data, it has to have the 'y(cm)', 'mean' and 'std. dev.' columns
            ylabel (str): string to write in the plot ylabel (ax.set_ylabel() will be used)
            dtype_label (str, optional): string printed in the plot for possible additional information (e.g. tally type or reaction type). Defaults to ''.
        """

        self.measured_data = measured_data
        self.tickers = np.arange(len(measured_data))

        # get all values in plot-friendly format
        measured_rstd = rel_std_dev(measured_data)

        # some data are already strings other have to be decoded
        try:
            my_xlabels = [el.decode() for el in measured_data['y(cm)']]
        except AttributeError:
            my_xlabels = measured_data['y(cm)']

        # get data floor and ceiling for plot framing
        floor, ceiling = _get_floor_ceiling([measured_data['mean']])

        # plot
        self.ax1.set_yscale('log')
        self.ax1.set_ylim([floor, ceiling])
        self.ax1.set_xticks(self.tickers)
        self.ax1.tick_params(axis='x', labelbottom=False)
        self.ax1.tick_params(axis='both', which='both', direction='in')
        self.ax1.set_ylabel(ylabel, fontsize=12)
        # add errorbars as grey areas and 1,2,3sigma strings
        _3sigma_area(self.ax2, self.tickers, measured_rstd)
        _print_3sigma(
            self.ax2, measured_data['mean'], measured_data['std. dev.'], self.tickers)
        # add C/E=1 line for showing target value
        self.ax2.hlines(1.0, -1, 20, colors='k', linestyles='-',
                        linewidth=1, label='_nolegend_')
        self.ax2.set_xlim([-0.5, self.tickers[-1] + .6])
        self.ax2.set_ylim([0.5, 1.5])
        self.ax2.set_xticks(self.tickers)
        self.ax2.set_xticklabels(my_xlabels)
        self.ax2.tick_params(axis='x', labelrotation=45)
        self.ax2.tick_params(axis='both', which='both', direction='in')
        self.ax2.set_xlabel('Position (cm)', fontsize=12)
        self.ax2.set_ylabel('C/E', fontsize=12)
        self.ax2.annotate(dtype_label, [0.02, 0.07], xycoords='axes fraction',
                          horizontalalignment='left', verticalalignment='bottom', fontsize=12)

        self.ax1.plot(self.tickers, self.measured_data['mean'], marker='s', ms=10,
                      ls='none', mew=1.5, mec='k', mfc='none', alpha=1, label='Experiment')

    def add_computed_data(self, dataset, marker='o', color='tab:red', label=''):
        """Includes computed data in the already existing plot generated in the class initialization

        Args:
            dataset (dataframe or dataset): table with the computed data, it has to have the 'y(cm)', 'mean' and 'std. dev.' columns
            marker (str, optional): matplotlib marker type. Defaults to 'o'.
            color (str, optional): matplotlib color. Defaults to 'tab:red'.
            label (str, optional): matplotlib label for legend. Defaults to ''.
        """

        # try/except statements are necessary in the case the user did not do every type of simulation but just some
        try:
            rstd = rel_std_dev(dataset)
            ce = mean_ratio(dataset['mean'], self.measured_data['mean'])

            floor, ceiling = _get_floor_ceiling(
                [self.measured_data['mean'], dataset['mean']])

            self.ax1.set_ylim([floor, ceiling])
            self.ax1.plot(self.tickers, dataset['mean'], marker=marker,
                          ms=7, ls='none', alpha=1, color=color, label=label)
            self.ax2.errorbar(self.tickers, ce, rstd*ce, marker=marker, ms=6, capsize=4,
                              barsabove=True, zorder=9, color=color, ls='none', label='_label')
            self.ax1.legend(frameon=True, fontsize=12)
        except:
            pass


def get_openmc_tally(simulation_type: str, statepoint_file: str, tally_name: str):
    """Read automatically every openmc tally from a give statepoint file for a given simulation type

    Args:
        simulation_type (str): can be 'onaxis', 'offaxis' or 'heating'
        statepoint_file (str): name of the statepoint file to read
        tally_name (str): string with the exact name for the tally defined in the openmc model

    Raises:
        ValueError: the name of the simulation_type is not one of the acceptable ones

    Returns:
        dataframe: dataframe of the read tally
        None: if there was no result from openmc for the given simulation type
    """

    # raise error if wrong simulation_type
    if simulation_type not in _sim_types:
        raise ValueError(
            "Invalid simulation_type. Expected one of: %s" % _sim_types)

    # enter the right folder according to the simulation type, read the openmc statepoint file and convert to dataframe
    try:
        os.chdir(simulation_type)
        sp = openmc.StatePoint(statepoint_file)
        tally = sp.get_tally(name=tally_name).get_pandas_dataframe()
        os.chdir('..')
        return tally
    except FileNotFoundError:
        pass


def postprocess_openmc_foils(tally_dataframe, foil_group: str):
    """Postprocesses the openmc results for the reaction rates simulation type (either onaxis and offaxis)

    Args:
        tally_dataframe (dataframe): dataframe with the tally results from openmc
        foil_group (str): can be "onaxis_group1', 'onaxis_group2' or 'offaxis'

    Raises:
        ValueError: the name of the foil_group is not one of the acceptable ones

    Returns:
        dataframe: dataframe with correctly postprocessed results from openmc (i.e. results normalized over correct foil volumes)
    """

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

    return df


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
    omc_qtld = [omc_qn[i] * cn * ce + qp * cp for i, qp in enumerate(omc_qp)]

    # rewrite the tally results dataframe now consisten with measured data
    tally_dataframe = tally_dataframe.drop(
        columns=['particle', 'nuclide']).drop_duplicates('cell')
    tally_dataframe['mean'] = omc_qtld[0]
    tally_dataframe['std. dev.'] = omc_qtld[1]

    return tally_dataframe
