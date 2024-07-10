import numpy as np
from typing import Iterable
import math
import matplotlib.axes
import matplotlib.pyplot as plt
from abc import ABC


def add_floor_ceiling(ax: matplotlib.axes, values: Iterable, scale: str = 'lin', gap: float = 0.):
    """This function computes the minimum and maximum values of a set of different arrays
    collected in a single list. It gets useful for finding the y_limits of a plot when all
    the values plotted are not known a priori

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
    def __init__(self, figsize, height_ratios, xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        self._xaxis = xaxis
        self._ylabel = ylabel
        self._dtype_label = dtype_label

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

    def add_reference_results(self, reference_data):
        self.reference_data = reference_data
        self._reference_tickers = np.arange(len(self.reference_data))

    def add_computed_results(self, computed_data):
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
        plt.show()


class PlotActivationFoils(PlotResults):
    def __init__(self, figsize=(6, 5), height_ratios=[2, 1.25], xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        super().__init__(figsize, height_ratios, xaxis, ylabel, dtype_label)

    def add_reference_results(self, reference_data, marker='s', color='k', alpha=1., label=''):
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
                          linewidth=1, label='_nolegend_')
        for ax in self.ax:
            ax.set_xlim([self._reference_tickers[0]-0.5,
                        self._reference_tickers[-1] + .6])
            ax.set_xticks(self._reference_tickers)
        self.ax[1].set_ylim([0.1, 1.75])
        self.ax[1].set_xticklabels(self.reference_data[self._xaxis])

        self.ax[0].legend(frameon=True, fontsize=12)

    def add_computed_results(self, computed_data, marker='o', color='tab:red', alpha=1., label=''):
        super().add_computed_results(computed_data)

        self.ax[0].plot(self._computed_tickers, computed_data['mean'], marker=marker,
                        ms=7, ls='none', alpha=alpha, color=color, label=label)
        self.ax[1].errorbar(self._computed_tickers, self.ce, self.rstd*self.ce, marker=marker, ms=6, capsize=4,
                            barsabove=True, zorder=9, color=color, ls='none', alpha=alpha, label='_')
        self.ax[0].legend(frameon=True, fontsize=12)


class PlotEnergySpectra:
    def __init__(self, figsize=(6, 5), height_ratios=[2, 1.25], xaxis: str = '', ylabel: str = '', dtype_label: str = ''):
        super().__init__(figsize, height_ratios, xaxis, ylabel, dtype_label)

        self.ax[1].set_yscale('log')

    def add_reference_results(self, reference_data):
        super().add_reference_results(reference_data)

    def add_computed_results(self, computed_data):
        super().add_computed_results(computed_data)
