import numpy as np
from typing import Iterable
import math.floor
import matplotlib.axes


def get_floor_ceiling(values: Iterable, scale: str = 'lin', gap: float = 0.):
    """This function computes the minimum and maximum values of a set of different arrays
    collected in a single list. It gets useful for finding the y_limits of a plot when all
    the values plotted are not known a priori.

    Args:
        values (Iterable): iterable of arrays/lists of floats. The function finds the absolute
        min and max values of these arrais.
        scale (str, optional): can be "lin" or "log". If "lin" is selected the function returns
        the absolute min and max values present in the values iterable.
        If "log" is selected the function returns 1eX, 1eY where X and Y are the orders of
        magnitude of the absolute min and max. Defaults to "lin".
        gap (float, optional): In order to better frame a plot the function can return values 
        that are sure to embed all the data in the values argument. If gap=X is selected, the 
        function returns min-X and max+X as values when scale='lin'. In the case
        scale='log' is selected the gap argument gets subtracted and added to the order of 
        magnitude of min and max respectively. The user can then frame a loglog or semilogy plot
        directly changing the orders of magnitude of the two ylims. defaults to 0. 

    Raises:
        NameError: if scale is neither "lin" or "log" raises a NameError.

    Returns:
        float, float: min-gap, max+gap if scale='lin'
        10**(X-gap), 10**(Y+gap) where X is min's order of magnitude and Y is max's order of magnitude
        if scale='log'
    """

    if scale not in ['lin', 'log']:
        msg = f"Wrong scale argument. It must be either 'lin' or 'log'"
        raise NameError(msg)

    min_value = min([np.nanmin(i[np.nonzero(i)]) for i in values])
    max_value = max([np.nanmax(i[np.nonzero(i)]) for i in values])

    if scale == 'lin':

        return min_value - gap, max_value + gap

    elif scale == 'log':
        min_oom = math.floor(math.log(min_value, 10))
        max_oom = math.floor(math.log(max_value, 10))

        return 10**(min_oom-gap),  10**(max_oom+gap)


def plot_stddev_area(ax: matplotlib.axes, ticks: Iterable, mean: Iterable, std_dev: Iterable,
                     color: str = 'k', alpha: float = .1, uncertainty_deg: int = 3):
    """This function the standard deviations of a set of data as shaded areas to a plot that has to
    already exist. it is possible to chose whether to plot 1, 2 or 3 times the standard deviations.
    It is based on the matplitlib.axes.fill_between() function.

    Args:
        ax (matplotlib.axes): the matplotlib.axes object of the plot that has already to exist
        ticks (Iterable): plot ticks for the x-axis
        mean (Iterable): list or array of the data statistical means
        std_dev (Iterable): list or array of the standard deviations to plot (in absolute value, not relative)
        color (str, optional): std. dev. areas' color, according to matplotlib. Defaults to 'k'.
        alpha (float, optional): std dev. areas' transparency degree according to matplotlib. Defaults to .1.
        uncertainty_deg (int, optional): Integer that can be 1, 2 or 3. Describes how many std. dev. to plot
        Defaults to 3.

    Raises:
        ValueError: if the uncertainty_deg arg is not 1, 2 or 3 raises a ValueError.
    """

    if uncertainty_deg not in [1, 2, 3]:
        msg = f'Value {uncertainty_deg} is not valid. It has to be an integer in [1, 2, 3]'
        raise ValueError(msg)

    ax.fill_between(ticks, mean - std_dev, mean +
                    std_dev, color=color, alpha=alpha)
    if uncertainty_deg > 1:
        ax.fill_between(ticks, mean - 2*std_dev, mean + 2 *
                        std_dev, color=color, alpha=alpha)
    if uncertainty_deg == 3:
        ax.fill_between(ticks, mean + 3*std_dev, mean - 3 *
                        std_dev, color=color, alpha=alpha)
