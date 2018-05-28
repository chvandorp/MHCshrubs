from __future__ import(print_function, division)
from builtins import (map, zip)
import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
from mhcshrubs import definitions
import mhcshrubs.auxiliary as aux

def CImeanBootstrap(xs, mass=0.95, B=1000):
    N = len(xs)
    ms = [np.mean(np.random.choice(xs, size=N)) for _ in range(B)]
    lo = np.percentile(ms, 100*(1.0-mass)/2)
    hi = np.percentile(ms, 100*(1.0 - (1.0-mass)/2))
    return (lo, hi)

def calcWAIC(loglikes, verbose=False, mass=0.95):
    """
    Compute the WAIC from the MCMC output.

    Args:
        loglikes -- for each MCMC sample, a vector of the log-likelihood
            of each observation

    Kwargs:
        verbose (bool) -- print the result, and some other info.
        mass (float) -- used for calculating the CI

    Returns:
        WAIC, WAIC_se, p_waic_hat, lpd_hat, WAIC_lo, WAIC_hi
        where [WIAC_lo, WAIC_hi] is the X% confidence range with X = 100 * mass
    """
    loglikes = [list(map(float, lls)) for lls in aux.transpose(loglikes)]
    lpd_hat_vec = [np.log(np.mean(list(map(np.exp, lls)))) for lls in loglikes]
    p_waic_hat_vec = [np.var(lls, ddof=1) for lls in loglikes]
    elpd_waic_hat_vec = [x-y for x, y in zip(lpd_hat_vec, p_waic_hat_vec)]
    WAIC_vec = np.array([-2*x for x in elpd_waic_hat_vec])
    WAIC = np.sum(WAIC_vec)
    N = len(WAIC_vec)
    WAIC_se = np.sqrt(N) * np.std(WAIC_vec)
    WAIC_lo, WAIC_hi = CImeanBootstrap(N * WAIC_vec, mass=mass)
    p_waic_hat = np.sum(p_waic_hat_vec)
    lpd_hat = np.sum(lpd_hat_vec)
    if verbose:
        print("number of observations = {}".format(N))
        print("WAIC = {:0.2f}".format(WAIC))
        print("standard error = {:0.2f}".format(WAIC_se))
        print("effective number of parameters = {:0.2f}".format(p_waic_hat))
        print("log pointwise predictive density = {:0.2f}".format(lpd_hat))
    return (WAIC, WAIC_se, p_waic_hat, lpd_hat, WAIC_lo, WAIC_hi)


def calcWBIC(loglike, verbose=False):
    """
    returns WBIC. NB: this is only correct when the sampling temperature
    was correctly set
    """
    WBIC = -2 * np.mean(loglike)
    print("WBIC = {:.2f}".format(WBIC))
    return WBIC


def calcGelmanRubinRhat(chains, parname):
    """
    Compute the Gelman-Rubin diagnostic for convergence of the MCMC

    Args:
        chains (list) -- a list of chains (a chain is a dict keyed by the parameter names)
        parname (str) -- the name of the parameter of interest

    Returns:
        Rhat -- the Gelman-Rubin statistic
    """
    traces = [chain[parname] for chain in chains]
    n = len(traces[0]) ## length of the chain(s)
    thetabars = np.array([np.mean(trace, axis=0) for trace in traces])
    variances = np.array([np.var(trace, axis=0, ddof=1) for trace in traces])
    W = np.mean(variances, axis=0)
    B = n * np.var(thetabars, axis=0, ddof=1)
    varhat = (1-1/n)*W + B/n
    with np.errstate(divide='ignore', invalid='ignore'):
        Rhat = np.sqrt(varhat/W)
    return Rhat

def mkTracePlot(ax, chains, parname, symbol=None, idx=0):
    """
    make a trace plot of a parameter, using multiple chains.

    Args:
        ax: pyplot object
        chains: a list of MCMCs
        parname: the name of the parameter of interest

    Kwargs:
        symbol: symbol used for label. If None (default), parname is used.
        idx: if the parameter is vector-valued, an index can be given.
            Otherwise, the first element of the vector is used

    Returns:
        None -- the ax object is modified
    """
    traces = [chain[parname] for chain in chains]
    ## if nothing is pased, assume scalar
    scalar = not isinstance(traces[0][0], list) if len(traces) > 0 else True
    for trace in traces:
        ts = range(len(trace))
        xs = [x if scalar else x[idx] for x in trace]
        ax.plot(ts, xs, alpha=np.sqrt(1.0/len(traces)))
    ax.set_xlabel("state")
    if scalar:
        if symbol is None:
            ax.set_ylabel(parname)
        else:
            ax.set_ylabel(symbol)
    else:
        if symbol is None:
            ax.set_ylabel("{0}[{1}]".format(parname, idx))
        else:
            if len(symbol) > 2 and symbol[0] == '$' and symbol[-1] == '$':
                ## this is TeX: remove the $s before adding the subscript
                ax.set_ylabel("${{{0}}}_{{{1}}}$".format(symbol[1:-1], idx))
            else:
                ax.set_ylabel("{0}[{1}]".format(symbol, idx))

def mkTracePlots(filename, chains, parnames, Rhat={}, cols=3):
    """
    make a number of trace plots and render a figure.

    Args:
        filename (str): the path to the figure
        chains: a list of MCMCs
        parnames: a list of pairs (parname (str), index (int or None))
            where indices is None indicates that the parameter is scalar-valued

    Kwargs:
        Rhat: a dictionary of parname -> Gelman-Rubin R hat diagnostic.
            parameters not in Rhat are ignored (default: empty dict)
        cols: the number of columns (default: 3)

    Returns:
        None -- a figure is created
    """
    cols = 3
    rem = len(parnames) % cols
    rows = len(parnames) // cols + (0 if rem == 0 else 1)
    fig, axs = plt.subplots(rows, cols, figsize=(15,3*rows))
    for i, pn in enumerate(parnames):
        col = i % cols
        row = i // cols
        symb = definitions.symbolDict[pn[0]] if pn[0] in definitions.symbolDict.keys() else None
        mkTracePlot(axs[row, col], chains, pn[0], symbol=symb, idx=pn[1])
        if pn[0] in Rhat.keys():
            R = Rhat[pn[0]][pn[1]] if pn[1] is not None else Rhat[pn[0]]
            axs[row, col].set_title("$\\hat{{R}} = {:0.2f}$".format(R))
    fig.tight_layout()
    fig.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig) ## prevent warnings when running in parallel
