import numpy as np
from typing import Iterable
import math
import matplotlib


def get_floor_ceiling(values: Iterable, scale: str = 'lin', gap: float = 0.):

    if scale not in ['lin', 'log']:
        msg = f"Wrong scale argument. It must be either 'lin' or 'log'"
        raise NameError(msg)

    min_value = np.nanmin([np.nanmin(i) for i in values])
    max_value = np.nanmax([np.nanmax(i) for i in values])

    if scale == 'lin':

        return min_value - gap, max_value + gap

    elif scale == 'log':
        min_oom = math.floor(math.log(min_value, 10))
        max_oom = math.floor(math.log(max_value, 10))

        return 10**(min_oom-gap),  10**(max_oom+gap)


def plot_stddev_area(ax: matplotlib.axes, ticks: Iterable, mean: Iterable, std_dev: Iterable, color: str = 'k', alpha: float = .1, uncertainty_deg: int = 3):

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
