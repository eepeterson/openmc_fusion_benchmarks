import numpy as np
from scipy.stats import norm


def get_gauss(mu, sigma):
    x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)
    y = norm.pdf(x, mu, sigma)
    y = y / max(y)

    return x, y
