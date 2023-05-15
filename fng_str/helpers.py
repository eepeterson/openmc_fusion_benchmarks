import numpy as np
import math
import matplotlib.pyplot as plt


def _get_floor_ceiling(values):
    # min value and max value for rr floor and ceiling

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
    ax.fill_between(ticks, 1 - std_dev, 1 + std_dev, color='k', alpha=.1)
    ax.fill_between(ticks, 1 - 2*std_dev, 1 + 2*std_dev, color='k', alpha=.1)
    ax.fill_between(ticks, 1 - 3*std_dev, 1 + 3*std_dev, color='k', alpha=.1)


def _print_3sigma(ax, mean, std_dev, ticks):
    # get index of last valid value for measured data for placing 1,2,3sigma strings near the grey areas
    sigma_index = (~np.isnan(mean)).cumsum(0).argmax(0)
    last_sigma = std_dev[sigma_index]

    if last_sigma < .08:
        last_sigma = .08
        sigma_start = .93

    ax.annotate('1\u03C3', [ticks[sigma_index] + 0.07,
                1. + 1*last_sigma], fontsize=10, clip_on=False)
    ax.annotate('2\u03C3', [ticks[sigma_index] + 0.07,
                1. + 2*last_sigma], fontsize=10, clip_on=False)
    ax.annotate('3\u03C3', [ticks[sigma_index] + 0.07,
                1. + 3*last_sigma], fontsize=10, clip_on=False)


class VisualizeResults:

    def __init__(self):

        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6, 5),
                                       gridspec_kw={'height_ratios': [2, 1.25]}, constrained_layout=True)

        self.fig = fig
        self.ax1 = ax1
        self.ax2 = ax2

    def add_measured_data(self, measured_data, tickers, ylabel: str, dtype_label=''):

        self.measured_data = measured_data

        # get all values in plot-friendly format
        measured_rstd = rel_std_dev(measured_data)

        my_xlabels = measured_data['y(cm)']
        # min value and max value for rr floor and ceiling
        floor, ceiling = _get_floor_ceiling([measured_data['mean']])

        # plot

        self.ax1.set_yscale('log')
        self.ax1.set_ylim([floor, ceiling])
        self.ax1.set_xticks(tickers)
        self.ax1.tick_params(axis='x', labelbottom=False)
        self.ax1.tick_params(axis='both', which='both', direction='in')
        self.ax1.set_ylabel(ylabel, fontsize=12)
        _3sigma_area(self.ax2, tickers, measured_rstd)
        _print_3sigma(
            self.ax2, measured_data['mean'], measured_data['std. dev.'], tickers)
        self.ax2.hlines(1.0, -1, 20, colors='k', linestyles='-',
                        linewidth=1, label='_nolegend_')
        self.ax2.set_xlim([-0.5, tickers[-1] + .6])
        self.ax2.set_ylim([0.5, 1.5])
        self.ax2.set_xticks(np.arange(len(tickers)))
        self.ax2.set_xticklabels(my_xlabels)
        self.ax2.tick_params(axis='x', labelrotation=45)
        self.ax2.tick_params(axis='both', which='both', direction='in')
        self.ax2.set_xlabel('Position (cm)', fontsize=12)
        self.ax2.set_ylabel('C/E', fontsize=12)
        self.ax2.annotate(dtype_label, [0.02, 0.07], xycoords='axes fraction',
                          horizontalalignment='left', verticalalignment='bottom', fontsize=12)

        self.ax1.plot(tickers, self.measured_data['mean'], marker='s', ms=10,
                      ls='none', mew=1.5, mec='k', mfc='none', alpha=1, label='Experiment')

    def add_computed_data(self, dataset, tickers, marker='o', color='tab:red', label=''):

        rstd = rel_std_dev(dataset)
        ce = mean_ratio(dataset['mean'], self.measured_data['mean'])

        floor, ceiling = _get_floor_ceiling(
            [self.measured_data['mean'], dataset['mean']])

        self.ax1.set_ylim([floor, ceiling])
        self.ax1.plot(tickers, dataset['mean'], marker=marker,
                      ms=7, ls='none', alpha=1, color=color, label=label)
        self.ax2.errorbar(tickers, ce, rstd*ce, marker=marker, ms=6, capsize=4,
                          barsabove=True, zorder=9, color=color, ls='none', label='_label')
        self.ax1.legend(frameon=True, fontsize=12)
