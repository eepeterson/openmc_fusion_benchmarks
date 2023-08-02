import numpy as np
import math
import matplotlib.pyplot as plt
import os

import openmc

# data necessary for postprocessing
# reaction rate simulations
# openmc model foil cell id per group
_foil_cells = [301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311]
# foil cell volumes per group
_detector_volume = 4.**2


# functions for postprocessing and visualization
def _get_floor_ceiling(values):
    """Finds the global minimum and maximum values that can be used to frame a plot

    Args:
        values (array): values that will be plotted

    Returns:
        float, float: global minimum and global maximum values present in the input array
    """

    min_value = min([np.nanmin(i[np.nonzero(i)]) for i in values])
    max_value = max([np.nanmax(i[np.nonzero(i)]) for i in values])
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


def lethargy(dataset):
    return np.log(dataset['energy high [eV]']/dataset['energy low [eV]'])


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


class PlotReactionRates:
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
            my_xlabels = [el.decode() for el in measured_data['pos. no.']]
        except AttributeError:
            my_xlabels = measured_data['pos. no.']

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
                [self.measured_data['mean'], np.array(dataset['mean'])])

            self.ax1.set_ylim([floor, ceiling])
            self.ax1.plot(self.tickers, dataset['mean'], marker=marker,
                          ms=7, ls='none', alpha=1, color=color, label=label)
            self.ax2.errorbar(self.tickers, ce, rstd*ce, marker=marker, ms=6, capsize=4,
                              barsabove=True, zorder=9, color=color, ls='none', label='_label')
            self.ax1.legend(frameon=True, fontsize=12)
        except:
            pass


class PlotEnergySpectra:
    """Embeds all necessary functions for plotting results and C/E results for comparison
    __init__() already generates the plot istance
    """

    def __init__(self):

        fig, ax = plt.subplots()

        self.fig = fig
        self.ax = ax

    def add_measured_data(self, measured_data, ylabel: str, dtype_label=''):

        x = measured_data['energy low [eV]']
        measured_mean = mean(measured_data) / lethargy(measured_data)
        measured_stddev = std_dev(measured_data) / lethargy(measured_data)

        self.ax.step(x, measured_mean, lw=1.5, c='k', label='Experiment')
        self.ax.fill_between(x, measured_mean - measured_stddev, measured_mean +
                             measured_stddev, step='pre', color='k', alpha=0.2)
        self.ax.set_yscale('log')
        self.ax.set_xlabel('Energy (eV)', fontsize=12)
        self.ax.set_ylabel(ylabel, fontsize=12)
        self.ax.tick_params(axis='both', which='both',
                            direction='in', labelsize=12)
        self.ax.legend(fontsize=12)

    def add_computed_data(self, dataset, color='tab:red', label=''):

        x = dataset['energy low [eV]']
        data_mean = mean(dataset) / lethargy(dataset)
        data_stddev = std_dev(dataset) / lethargy(dataset)

        try:
            self.ax.step(x, data_mean, lw=1.5, c=color, label=label)
            self.ax.fill_between(x, data_mean - data_stddev, data_mean +
                                 data_stddev, step='pre', color=color, alpha=0.2)
            self.ax.legend(frameon=True, fontsize=12)
        except:
            pass


def get_openmc_tally(statepoint_file: str, tally_name: str):
    """Read automatically every openmc tally from a give statepoint file for a given simulation type

    Args:
        statepoint_file (str): name of the statepoint file to read
        tally_name (str): string with the exact name for the tally defined in the openmc model

    Raises:
        ValueError: the name of the simulation_type is not one of the acceptable ones

    Returns:
        dataframe: dataframe of the read tally
        None: if there was no result from openmc for the given simulation type
    """

    # enter the right folder according to the simulation type, read the openmc statepoint file and convert to dataframe
    # try:
    os.chdir('results')
    sp = openmc.StatePoint(statepoint_file)
    tally = sp.get_tally(name=tally_name).get_pandas_dataframe()
    os.chdir('..')
    return tally
    # except FileNotFoundError:
    #     pass


def postprocess_openmc_foils(tally_dataframe):
    """Postprocesses the openmc results for the reaction rates simulation type

    Args:
        tally_dataframe (dataframe): dataframe with the tally results from openmc

    Returns:
        dataframe: dataframe with correctly postprocessed results from openmc (i.e. results normalized over correct foil volumes)
    """

    df = tally_dataframe.loc[(tally_dataframe['particle'] == 'neutron')]

    df = df[df['cell'].isin(_foil_cells)]
    df['mean'] = df['mean'] / _detector_volume
    df['std. dev.'] = df['std. dev.'] / _detector_volume

    return df


def postprocess_openmc_spectrum(tally_dataframe, cell):
    """Postprocesses the openmc results for the neutron energy spectra

    Args:
        tally_dataframe (dataframe): dataframe with correctly postprocessed results from openmc (i.e. results normalized over lethargy and detector cell volumes)
    """

    df = tally_dataframe.loc[(
        tally_dataframe['particle'] == 'neutron')]

    df = df.loc[(df['cell'] == cell)]
    df['mean'] = df['mean'] / _detector_volume
    df['std. dev.'] = df['std. dev.'] / _detector_volume

    return df
