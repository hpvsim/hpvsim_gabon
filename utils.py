"""
Utilities
"""

# Imports
import sciris as sc
import numpy as np
from scipy.stats import norm, lognorm


def set_font(size=None, font='Libertinus Sans'):
    """ Set a custom font """
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


def logn_percentiles_to_pars(x1, p1, x2, p2):
    """ Find the parameters of a lognormal distribution where:
            P(X < p1) = x1
            P(X < p2) = x2
    """
    x1 = np.log(x1)
    x2 = np.log(x2)
    p1ppf = norm.ppf(p1)
    p2ppf = norm.ppf(p2)
    s = (x2 - x1) / (p2ppf - p1ppf)
    mean = ((x1 * p2ppf) - (x2 * p1ppf)) / (p2ppf - p1ppf)
    scale = np.exp(mean)
    return s, scale


def get_debut(sex='f'):
    """
    Read in dataframes taken from DHS and return them in a plot-friendly format,
    optionally saving the distribution parameters

    Percentiles below are derived by fitting to the 2019-21 Gabon DHS; see
    https://www.researchsquare.com/article/rs-3074559/v1
    """
    if sex == 'f':
        x1 = 15
        p1 = 0.184
        x2 = 20
        p2 = 0.858

    else:
        x1 = 15
        p1 = 0.203
        x2 = 20
        p2 = 0.786
    s, scale = logn_percentiles_to_pars(x1, p1, x2, p2)
    rv = lognorm(s=s, scale=scale)

    return rv.mean(), rv.std()


def plot_single(ax, mres, to_plot, start_year, end_year, color, ls='-', label=None, smooth=True, smooth_window=5):
    """Plot one metric's median + envelope over [start_year, end_year], optionally smoothed.

    ``mres`` is ``msim_dict[scen]`` — a dict-like keyed by metric name, where each
    value is a pandas DataFrame indexed by year with columns ``median``, ``low``,
    ``high`` (built by ``run_scenarios.run_sims``).
    """
    df = mres[to_plot]
    sub = df.loc[(df.index >= start_year) & (df.index <= end_year)]
    years = sub.index.to_numpy()
    best = sub['median'].to_numpy()
    low = sub['low'].to_numpy()
    high = sub['high'].to_numpy()

    if smooth:
        best = np.convolve(best, np.ones(smooth_window), 'valid') / smooth_window
        low = np.convolve(low, np.ones(smooth_window), 'valid') / smooth_window
        high = np.convolve(high, np.ones(smooth_window), 'valid') / smooth_window
        years = years[smooth_window - 1:]

    ax.plot(years, best, color=color, label=label, ls=ls)

    if to_plot == 'asr_cancer_incidence':
        elim_year = sc.findfirst(best < 4, die=False)
        if elim_year is not None:
            print(f'{label} elim year: {years[elim_year]}')
        else:
            print(f'{label} not eliminated')

    ax.fill_between(years, low, high, alpha=0.1, color=color)

    # WHO elimination threshold (4 per 100,000)
    ax.axhline(4, color='k', ls='--', lw=0.5)
    return ax
    return ax


# %% Run as a script
if __name__ == '__main__':

    for sex in ['f', 'm']:
        mean, std = get_debut(sex)
        print(f'Mean debut age ({sex}): {mean:.2f}, std: {std:.2f}')

