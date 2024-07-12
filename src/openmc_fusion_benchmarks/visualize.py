import numpy as np
import pandas as pd
from typing import Iterable
import math
import matplotlib.axes
import matplotlib.pyplot as plt
from abc import ABC

import openmc_fusion_benchmarks as ofb


def add_floor_ceiling(ax: matplotlib.axes, values: Iterable, scale: str = 'lin', gap: float = 0.):
    """This function computes the minimum and maximum values of a set of different arrays
    collected in a single list. It gets useful for finding the y_limits of a plot when all
    the values plotted are not known a priori.

    Parameters
    ----------
    values : Iterable
        iterable of arrays/lists of floats. The function finds the absolute
        min and max values of these arrays
    scale : str, optional
        can be "lin" or "log". If "lin" is selected the function returns the
        absolute min and max values present in the values iterable.
        If "log" is selected the function returns 1eX, 1eY where X and Y are the orders of
        magnitude of the absolute min and max, by default 'lin'
    gap : float, optional
        In order to better frame a plot the function can return values
        that are sure to embed all the data in the values argument. If gap=X is selected, the
        function returns min-X and max+X as values when scale='lin'. In the case
        scale='log' is selected the gap argument gets subtracted and added to the order of
        magnitude of min and max respectively. The user can then frame a loglog or semilogy plot
        directly changing the orders of magnitude of the two ylims, by default 0.

    Returns
    -------
    float, float
        min-gap, max+gap if scale='lin' 10**(X-gap), 10**(Y+gap) where X is min's order of
        magnitude and Y is max's order of magnitude if scale='log'

    Raises
    ------
    NameError
        If scale is neither "lin" or "log" raises a NameError
    """
    # check scale argument is right
    if scale not in ['lin', 'log']:
        msg = f"Wrong scale argument. It must be either 'lin' or 'log'"
        raise NameError(msg)

    # get global min ang max values
    values = np.array(values)
    if isinstance(values, list):
        min_value = min([np.nanmin(i[np.nonzero(i)]) for i in values])
        max_value = max([np.nanmax(i[np.nonzero(i)]) for i in values])
    else:
        min_value = np.nanmin(values[np.nonzero(values)])
        max_value = np.nanmax(values[np.nonzero(values)])

    # return ylim([min, max]) in either linear or logarithmic form
    if scale == 'lin':

        ax.set_ylim(min_value - gap, max_value + gap)

    elif scale == 'log':
        min_oom = math.floor(math.log(min_value, 10))
        max_oom = math.floor(math.log(max_value, 10))

        ax.set_ylim(10**(min_oom-gap),  10**(max_oom+1+gap))


def plot_stddev_area(ax: matplotlib.axes, ticks: Iterable, mean: Iterable, std_dev: Iterable,
                     color: str = 'k', alpha: float = .1, uncertainty_deg: int = 3):
    """This function the standard deviations of a set of data as shaded areas to a plot that has to
    already exist. it is possible to chose whether to plot 1, 2 or 3 times the standard deviations.
    It is based on the matplitlib.axes.fill_between() function.

    Parameters
    ----------
    ax : matplotlib.axes
        the matplotlib.axes object of the plot that has already to exist
    ticks : Iterable
        plot ticks for the x-axis
    mean : Iterable
        list or array of the data statistical means
    std_dev : Iterable
        list or array of the standard deviations to plot (in absolute value, not relative)
    color : str, optional
        std. dev. areas' color, according to matplotlib, by default 'k'
    alpha : float, optional
        std dev. areas' transparency degree according to matplotlib, by default .1
    uncertainty_deg : int, optional
        Integer that can be 1, 2 or 3. Describes how many std. dev. to plot, by default 3

    Raises
    ------
    ValueError
        if the uncertainty_deg arg is not 1, 2 or 3 raises a ValueError
    """
    # raise error if uncertainty required is not between 1 and 3 sigma
    if uncertainty_deg not in [1, 2, 3]:
        msg = f'Value {uncertainty_deg} is not valid. It has to be an integer in [1, 2, 3]'
        raise ValueError(msg)

    # fill between for generating the shaded areas
    ax.fill_between(ticks, mean - std_dev, mean +
                    std_dev, color=color, alpha=alpha)
    if uncertainty_deg > 1:
        ax.fill_between(ticks, mean - 2*std_dev, mean + 2 *
                        std_dev, color=color, alpha=alpha)
    if uncertainty_deg == 3:
        ax.fill_between(ticks, mean + 3*std_dev, mean - 3 *
                        std_dev, color=color, alpha=alpha)


class PlotResults(ABC):
    """Abstract class for plotting results of the benchmarks. 
    It is meant to be inherited by other plotting classes specific for
    different tally types. It contains the basic methods to plot the results.
    """

    def __init__(self, figsize: Iterable[float, float], height_ratios: Iterable[float, float],
                 xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        """Constructor of the PlotResults class that initializes the figure
        and axes objects. It also sets the basic parameters for the plot.
        As it is meant to compare results (most likely computed results 
        agaist experimental results) it most likely has two rows of subfigures
        The top subfigure(s) is to compare absolute values of results while
        the bottom subfigure(s) is to show the C/E ratio of the results.

        Parameters
        ----------
        figsize : Iterable[float, float]
            Size of the figure to plot (width, height)
        height_ratios : Iterable[float, float]
            Ratios of the heights of the subfigures in the two rows mentioned
            above
        xaxis : str, optional
            Name of the result dataframe column meant to be the values of the 
            x-axis of the plot, by default ''
        ylabel : str, optional
            Name for the plot y-label, by default ''
        dtype_label : str, optional
            Overall name to identify the plot, by default ''
        """
        self._figsize = figsize
        self._height_ratios = height_ratios
        self._xaxis = xaxis
        self._ylabel = ylabel
        self._dtype_label = dtype_label

    def add_reference_results(self, reference_data: pd.DataFrame):
        """Method to add the reference results to the plot. Most likely the
        experimental results. It also sets the reference tickers for the plot
        and it is used as the "E" values in the "C/E" ratio plot. Each plot
        needs to have one and only one reference results set.

        Parameters
        ----------
        reference_data : pd.DataFrame
            DataFrame containing the reference results to plot
        """
        self.reference_data = reference_data
        self._reference_tickers = np.arange(len(self.reference_data))

    def add_computed_results(self, computed_data: pd.DataFrame):
        """Method to add the computed results to the plot. Most likely the
        results of the simulations. It also sets the computed tickers for the
        computed plot. Each plot can have as many computed results as needed."""

        self.computed_data = computed_data
        self._computed_tickers = np.arange(len(self.computed_data))

        try:
            self.reference_data
        except AttributeError:
            msg = """You have to add the reference data with the
            add_reference_results method before adding the computed data."""
            raise AttributeError(msg)

        self.rstd = computed_data['std. dev.']/computed_data['mean']
        self.ce = computed_data['mean']/self.reference_data['mean']

    def show(self):
        """Shows the plot. Just a plt.show() wrapper.
        """
        plt.show()

    def savefig(self, path: str, dpi: int = 300):
        """Saves the plot to a file. Just a plt.savefig() wrapper.

        Parameters
        ----------
        path : str
            Path to save the plot, including the filename and extension
        dpi : int, optional
            Quality of the image to save, by default 300
        """
        self.fig.savefig(path, dpi=dpi, bbox_inches="tight")


class PlotReactionRates(PlotResults):
    """Class for plotting the reaction rates results. Inherits from 
    PlotResults.
    """

    def __init__(self, figsize=(6, 5), height_ratios=[2, 1.25], xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        super().__init__(figsize, height_ratios, xaxis, ylabel, dtype_label)

        self.fig, self.ax = plt.subplots(nrows=2, ncols=1, figsize=figsize,
                                         gridspec_kw={'height_ratios': height_ratios}, constrained_layout=True)
        self.ax[0].set_yscale('log')
        self.ax[0].set_ylabel(self._ylabel, fontsize=12)
        self.ax[0].tick_params(axis='x', labelbottom=False)
        self.ax[0].tick_params(axis='both', which='both', direction='in')

        self.ax[1].annotate(self._dtype_label, [0.02, 0.07], xycoords='axes fraction',
                            horizontalalignment='left', verticalalignment='bottom', fontsize=12)
        self.ax[1].set_xlabel(self._xaxis, fontsize=12)
        self.ax[1].set_ylabel('C/E', fontsize=12)
        self.ax[1].tick_params(axis='x', labelrotation=45)
        self.ax[1].tick_params(axis='both', which='both', direction='in')

    def add_reference_results(self, reference_data, marker: str = 's', color: str = 'k', alpha: float = 1., label=''):
        """Method to add the reference results to a reaction rate plot.
        Alongside the reference_data argument, inherited by the PlotResults
        class, this method also takes the marker, color, alpha and label.

        Parameters
        ----------
        marker : str, optional
            Marker shape like in matplotlib, by default 's'
        color : str, optional
            Marker color like in matplotlib, by default 'k'
        alpha : float, optional
            Marker alpha of marker like in matplotlib, by default 1.
        label : str, optional
            Label for legend, by default ''
        """
        super().add_reference_results(reference_data)

        # plot results
        self.ax[0].plot(self._reference_tickers, self.reference_data['mean'], marker=marker, ms=10,
                        ls='none', mew=1.5, mec=color, mfc='none', alpha=alpha, label=label)

        plot_stddev_area(ax=self.ax[1], ticks=self._reference_tickers, mean=np.ones(len(
            self.reference_data['mean'])), std_dev=self.reference_data['std. dev.']/self.reference_data['mean'])
        add_floor_ceiling(ax=self.ax[0], values=[
                          self.reference_data['mean']], scale='log', gap=0.)
        #
        self.ax[1].hlines(1.0, self._reference_tickers[0]-1, self._reference_tickers[-1] + 1, colors='k', linestyles='-',
                          linewidth=1, label='_')
        for ax in self.ax:
            ax.set_xlim([self._reference_tickers[0]-0.5,
                        self._reference_tickers[-1] + .6])
            ax.set_xticks(self._reference_tickers)
        self.ax[1].set_ylim([0.1, 1.75])
        self.ax[1].set_xticklabels(self.reference_data[self._xaxis])

        self.ax[0].legend(frameon=True, fontsize=12)

    def add_computed_results(self, computed_data, marker: str = 'o', color: str = 'tab:red', alpha: float = 1., label=''):
        """Method to add the computed results to a reaction rate plot.
        Alongside the reference_data argument, inherited by the PlotResults
        class, this method also takes the marker, color, alpha and label.
        Which are described in the add_reference_results method.
        """
        super().add_computed_results(computed_data)

        self.ax[0].plot(self._computed_tickers, computed_data['mean'], marker=marker,
                        ms=7, ls='none', alpha=alpha, color=color, label=label)
        self.ax[1].errorbar(self._computed_tickers, self.ce, self.rstd*self.ce, marker=marker, ms=6, capsize=4,
                            barsabove=True, zorder=9, color=color, ls='none', alpha=alpha, label='_')
        self.ax[0].legend(frameon=True, fontsize=12)


class PlotCellHeating(PlotReactionRates):
    """Class for plotting the cell heating results. Inherits from 
    PlotReactionRates.
    """

    def __init__(self, figsize=(6, 5), height_ratios=[2, 1.25], xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        super().__init__(figsize, height_ratios, xaxis, ylabel, dtype_label)


class PlotEnergySpectra(PlotResults):
    """Class for plotting the energy spectra results. Inherits from
    PlotResults.
    """

    def __init__(self, figsize=(12, 6), height_ratios=[2, 1.25], xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        super().__init__(figsize, height_ratios, xaxis, ylabel, dtype_label)

        self.fig, self.ax = plt.subplots(nrows=2, ncols=2, figsize=figsize,
                                         gridspec_kw={'height_ratios': height_ratios}, constrained_layout=True)

        for i in range(2):
            self.ax[i, 0].set_xscale('log')
            self.ax[0, i].set_yscale('log')
            self.ax[1, i].set_yscale('log')
            self.ax[1, i].set_ylim([1e-2, 1e2])
            self.ax[0, i].tick_params(axis='x', labelbottom=False)
            self.ax[1, i].set_xlabel('Energy (eV)', fontsize=12)
            self.ax[1, i].xaxis.set_label_coords(0.5, -0.2)
            for j in range(2):
                self.ax[i, j].tick_params(
                    axis='both', which='both', direction='in', labelsize=12)

        self.ax[0, 0].set_ylabel(self._ylabel, fontsize=12)
        self.ax[1, 0].set_ylabel('C/E', fontsize=12)

        self.ax[1, 0].annotate("Log scale on x", [0.02, 0.07], xycoords='axes fraction',
                               horizontalalignment='left', verticalalignment='bottom', fontsize=12)
        self.ax[1, 1].annotate("Lin scale on x", [0.02, 0.07], xycoords='axes fraction',
                               horizontalalignment='left', verticalalignment='bottom', fontsize=12)

    def add_reference_results(self, reference_data, ls: str = '-', color: str = 'k', alpha: float = 1., label=''):
        """Method to add the reference results to an energy spectra plot.
        Alongside the reference_data argument, inherited by the PlotResults
        class, this method also takes the ls, color, alpha and label.

        Parameters
        ----------
        ls : str, optional
            Linestyle like in matplotlib, by default '-'
        color : str, optional
            Line color like in matplotlib, by default 'k'
        alpha : float, optional
            Line alpha like in matplotlib, by default 1.
        label : str, optional
            Label for legend, by default ''
        """
        super().add_reference_results(reference_data)

        min_ebound, max_ebound = ofb.get_nonzero_energy_interval(
            reference_data)
        min_oom = math.floor(math.log(min_ebound, 10))
        max_oom = math.floor(math.log(max_ebound, 10))

        for i in range(2):
            self.ax[0, i].step(reference_data['energy low [eV]'],
                               reference_data['mean'], ls=ls, lw=1.5, c=color, alpha=alpha, label=label)
            self.ax[0, i].fill_between(reference_data['energy low [eV]'], reference_data['mean'] - reference_data['std. dev.'], reference_data['mean'] +
                                       reference_data['std. dev.'], step='pre', color='k', alpha=.2*alpha)

            plot_stddev_area(ax=self.ax[1, i], ticks=reference_data['energy high [eV]'], mean=np.ones(len(
                reference_data['mean'])), std_dev=reference_data['std. dev.']/reference_data['mean'])

            self.ax[1, i].hlines(1.0, 0, np.array(reference_data['energy high [eV]'])[
                                 -1] + 5e6, colors='k', linestyles='-', linewidth=1, label='_')

            self.ax[i, 1].set_xlim([0, max_ebound + .1**max_oom])
            self.ax[i, 0].set_xlim(
                [max(min_ebound - .1*10**min_oom, .001), max_ebound + .1**max_oom])

        self.ax[0, 0].legend(frameon=True, fontsize=12)

    def add_computed_results(self, computed_data, ls='-', color='tab:red', alpha=1., label=''):
        """Method to add the computed results to an energy spectra plot.
        Alongside the computed_data argument, inherited by the PlotResults
        class, this method also takes the ls, color, alpha and label.
        Which are described in the add_reference_results method.
        """
        super().add_computed_results(computed_data)

        for i in range(2):
            self.ax[0, i].step(computed_data['energy low [eV]'],
                               computed_data['mean'], ls=ls, lw=1.5, c=color, alpha=alpha, label=label)
            self.ax[0, i].fill_between(computed_data['energy low [eV]'], computed_data['mean'] - computed_data['std. dev.'], computed_data['mean'] +
                                       computed_data['std. dev.'], step='pre', color=color, alpha=.2*alpha)

            self.ax[1, i].step(computed_data['energy low [eV]'],
                               self.ce, lw=1.5, c=color, alpha=alpha, label='_')

        self.ax[0, 0].legend(frameon=True, fontsize=12)
