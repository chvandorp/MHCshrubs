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

def mkAlleleWeightPlots(filename, chain, hlaAlleles):
    loci = sorted(list(hlaAlleles.keys()))
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    for loc, ax in zip(loci, axs):
        mkAlleleViolinPlot(ax, chain["beta{}".format(loc)], hlaAlleles[loc],
                           defn.locusColorDict[loc], sortfun=np.mean)
        ax.set_ylabel("weight HLA-{}".format(loc))
        ax.axhline(y=0, color='k')
        ax.set_ylim(-1,1)
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")

def mkAlleleFreqPlots(filename, chain, hlaAlleles, hlaCounts):
    loci = sorted(list(hlaAlleles.keys()))
    fig, axs = plt.subplots(len(loci), 1, figsize=(15, 3*len(loci)))
    for loc, ax in zip(loci, axs):
        mkAlleleViolinPlot(ax, chain["p{}".format(loc)], hlaAlleles[loc],
                           defn.locusColorDict[loc], alpha=0.5)
        freqs = [c/np.sum(hlaCounts[loc]) for c in hlaCounts[loc]] ## float division
        ax.scatter(range(len(freqs)), freqs, s=10, color='k', marker='D')
        ax.set_ylabel("frequency HLA-{}".format(loc))
    axs[-1].set_xlabel("allele")
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
