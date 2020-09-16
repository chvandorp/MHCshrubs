import numpy as np
import matplotlib.pyplot as plt

from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn

def mkAlleleViolinPlot(ax, trace, alleles, facecolor, edgecolor=None,
                       sortfun=None, alpha=1.0):
    """
    Make violin plots of the MCMC results
    """
    data = aux.transpose(trace)
    pos = range(len(data))
    sort_idxs = np.argsort([sortfun(x) for x in data]) if sortfun is not None else pos
    vp = ax.violinplot([data[i] for i in sort_idxs], pos,
                       showmeans=False, showextrema=False, showmedians=False)
    for pc in vp['bodies']:
        pc.set_facecolor(facecolor)
        pc.set_edgecolor(edgecolor)
        pc.set_alpha(alpha)
    ax.set_xticks(pos)
    ax.set_xticklabels([alleles[i].short_str() for i in sort_idxs], rotation=90)
    ax.set_xlim(-1, len(pos))

def mkAlleleWeightPlotsABC(filename, chain, hlaAlleles):
    loci = sorted(list(hlaAlleles.keys()))
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    if len(loci) == 1:
        axs = np.array([axs])
    for loc, ax in zip(loci, axs):
        mkAlleleViolinPlot(ax, chain["beta{}".format(loc)], hlaAlleles[loc],
                           defn.locusColorDict[loc], sortfun=np.mean)
        ax.set_ylabel("weight HLA-{}".format(loc))
        ax.axhline(y=0, color='k')
        ax.set_ylim(-1,1)
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")

def mkAlleleWeightPlots(filename, chain, hlaAlleles):
    loci = sorted(list(hlaAlleles.keys()))
    allelesPerLocus = [len(hlaAlleles[X]) for X in loci]
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    if len(loci) == 1:
        axs = np.array([axs])
    for ell, loc in enumerate(loci):
        ax = axs[ell]
        lidx, ridx = sum(allelesPerLocus[:ell]), sum(allelesPerLocus[:ell+1])
        ## the weights are concatenated in one array.
        betas = np.array(chain["beta"])[:,lidx:ridx].tolist()
        mkAlleleViolinPlot(ax, betas, hlaAlleles[loc],
                           defn.locusColorDict[loc], sortfun=np.mean)
        ax.set_ylabel("weight HLA-{}".format(loc))
        ax.axhline(y=0, color='k')
        ax.set_ylim(-1,1)
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")

def mkAlleleFreqPlotsABC(filename, chain, hlaAlleles, hlaCounts):
    loci = sorted(list(hlaAlleles.keys()))
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    if len(loci) == 1:
        axs = np.array([axs])
    for loc, ax in zip(loci, axs):
        mkAlleleViolinPlot(ax, chain["p{}".format(loc)], hlaAlleles[loc],
                           defn.locusColorDict[loc], alpha=0.5)
        sum_hla_counts = np.sum(hlaCounts[loc])
        if sum_hla_counts > 0: ## only plot HLA freqs if we have data...
            freqs = [c / sum_hla_counts for c in hlaCounts[loc]] ## x / y is float division
            ax.scatter(range(len(freqs)), freqs, s=10, color='k', marker='D')
        ax.set_ylabel("frequency HLA-{}".format(loc))
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")

def mkAlleleFreqPlots(filename, chain, hlaAlleles, hlaCounts):
    loci = sorted(list(hlaAlleles.keys()))
    allelesPerLocus = [len(hlaAlleles[X]) for X in loci]
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    if len(loci) == 1:
        axs = np.array([axs])
    for ell, loc in enumerate(loci):
        ax = axs[ell]
        lidx, ridx = sum(allelesPerLocus[:ell]), sum(allelesPerLocus[:ell+1])
        ps = np.array(chain["p"])[:,lidx:ridx].tolist()
        mkAlleleViolinPlot(ax, ps, hlaAlleles[loc],
                           defn.locusColorDict[loc], alpha=0.5)
        sum_hla_counts = np.sum(hlaCounts[loc])
        if sum_hla_counts > 0: ## only plot HLA freqs if we have data...
            freqs = [c / sum_hla_counts for c in hlaCounts[loc]] ## x / y is float division
            ax.scatter(range(len(freqs)), freqs, s=10, color='k', marker='D')
        ax.set_ylabel("frequency HLA-{}".format(loc))
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
