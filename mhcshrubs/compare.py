"""@package compare
Module for comparing models.
"""
from __future__ import (print_function, division)
from builtins import zip
import numpy as np
import scipy.stats as sts
import matplotlib.pyplot as plt
import csv
import os
from mhcshrubs import mhctools
from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn

def compareModels(modelname1, modelname2, deltalpd_threshold=0.1):
    """
    Compare two cross-validated models.
    """
    modelnames = [modelname1, modelname2]
    filenames = ["cv_alleles_summary.{}.tsv".format(modelname) for modelname in modelnames]
    work_folder = os.getcwd()
    fullfilenames = [os.path.join(work_folder, "data", filename) for filename in filenames]
    lpdDicts = [{}, {}]
    for filename, lpdDict in zip(fullfilenames, lpdDicts):
        with open(filename, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for record in reader:
                allele = mhctools.MhcObject(record["allele"])
                lpd = float(record["lpd"]) ## NB: using weighted lpd
                lpdDict[allele] = lpd
    alleles = sorted(lpdDicts[0].keys())
    loci = aux.unique([allele.locus for allele in alleles])
    alleles = dict((X, [allele for allele in alleles if allele.locus == X]) for X in loci)
    fig, axs = plt.subplots(1, len(loci), figsize=(15,10))
    height = np.max([len(alleles[locus]) for locus in loci])
    maxabsdiff = 0 ## used for x-limits of the axes
    for locus, ax in zip(loci, axs):
        selected_diffs = []
        selected_alleles = []
        for allele in alleles[locus]:
            vals = [lpdDict[allele] for lpdDict in lpdDicts]
            diff = vals[1] - vals[0]
            if np.abs(diff) >= deltalpd_threshold:
                selected_diffs.append(diff)
                selected_alleles.append(allele)
            maxabsdiff = max(np.abs(diff), maxabsdiff)
        pos = np.argsort(np.abs(selected_diffs))[::-1]
        ax.grid(axis='y')
        ax.barh(range(len(pos)), [selected_diffs[p] for p in pos], 0.3,
                color=defn.locusColorDict[locus])
        ax.set_yticks(range(len(pos)))
        ax.set_yticklabels([selected_alleles[p].short_str() for p in pos])
        ax.axvline(x=0, color='k')
        ax.set_xlabel("$\\Delta lpd$ HLA-{}".format(locus))
        ## do some basic statistics
        positive = np.sum(1 for x in selected_diffs if x > 0)
        total = len(selected_diffs)
        pval = sts.binom_test(positive, total)
        print("HLA-{0}: positive = {1}, total = {2}, P = {3}".format(locus, positive, total, pval))
    ## additional formatting...
    for ax in axs:
        ax.set_xlim((-maxabsdiff, maxabsdiff))
    fig.tight_layout()
    fig.savefig("compare.{0}.{1}.png".format(modelname1, modelname2), dpi=300, bbox_inches='tight')


def compareModelsSEM(modelname1, modelname2, deltalpd_threshold=0.1):
    """
    Compare two cross-validated models, also calculate the (scaled) SEMs
    """
    modelnames = [modelname1, modelname2]
    filenames = ["cv_alleles_summary.{}.tsv".format(modelname) for modelname in modelnames]
    work_folder = os.getcwd()
    fullfilenames = [os.path.join(work_folder, "data", filename) for filename in filenames]
    subjectDicts = [{}, {}]
    for filename, subjDict in zip(fullfilenames, subjectDicts):
        with open(filename, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for record in reader:
                allele = mhctools.MhcObject(record["allele"])
                subject_line = record["subjects"]
                subject_table = [row.split(",") for row in subject_line.split(";")]
                subject_dict = {int(row[0]) : (float(row[1]), float(row[2])) for row in subject_table}
                subjDict[allele] = subject_dict
    alleles = sorted(subjectDicts[0].keys())
    ## compute Delta lpd and scaled sems
    deltaLpdDict = {}
    semDict = {}
    for allele in alleles:
        subject_dicts = [subjDict[allele] for subjDict in subjectDicts]
        idxs = list(subject_dicts[0].keys())
        diffs = []
        weights = []
        for idx in idxs:
            clpd0, ppc0 = subject_dicts[0][idx]
            clpd1, ppc1 = subject_dicts[1][idx]
            if ppc0 != 0 and ppc1 != 0:
                diffs.append(clpd1 - clpd0)
                weights.append(ppc1 * ppc0)
        if len(diffs) > 0:
            deltaLpdDict[allele] = np.sum([d*w for d, w in zip(diffs, weights)])
            if len(diffs) > 1:
                mu = deltaLpdDict[allele] / np.sum(weights)
                semDict[allele] = np.sqrt(np.sum([w*(d-mu)**2 for d, w in zip(diffs, weights)]))
            else:
                semDict[allele] = np.nan
        else:
            deltaLpdDict[allele] = np.nan
            semDict[allele] = np.nan
    ## make figure
    loci = aux.unique([allele.locus for allele in alleles])
    alleles = dict((X, [allele for allele in alleles if allele.locus == X]) for X in loci)
    fig, axs = plt.subplots(1, len(loci), figsize=(10,6.67))
    height = np.max([len(alleles[locus]) for locus in loci])
    for locus, ax in zip(loci, axs):
        maxabsdiff = 0 ## used for x-limits of the axes
        selected_diffs = []
        selected_alleles = []
        selected_errs = []
        for allele in alleles[locus]:
            diff = deltaLpdDict[allele]
            if diff is not None and np.abs(diff) >= deltalpd_threshold:
                selected_diffs.append(diff)
                selected_alleles.append(allele)
                selected_errs.append(semDict[allele])
                maxabsdiff = max(np.abs(diff), maxabsdiff)
        ## revert order to get a "christmas tree plot"
        pos = np.argsort(np.abs(selected_diffs))[::-1]
        ax.grid(axis='y')
        ax.barh(range(len(pos)), [selected_diffs[p] for p in pos], 0.3,
                color=defn.locusColorDict[locus],
                xerr=np.array([selected_errs[p] for p in pos]),
                #ecolor=defn.locusColorDict[locus],
                capsize=5)
        ax.set_yticks(range(len(pos)))
        ax.set_yticklabels([selected_alleles[p].short_str() for p in pos])
        ax.axvline(x=0, color='k')
        ax.set_xlabel("$\\Delta {{\\rm lpd}}$ HLA-{}".format(locus))
        ## do some basic statistics
        positive = np.sum(1 for x in selected_diffs if x > 0)
        total = len(selected_diffs)
        pval_binom = sts.binom_test(positive, total)
        print("HLA-{0}: positive = {1}, total = {2}, P = {3} (binomial)".format(locus, positive, total, pval_binom))
        W, pval_wilcox = sts.wilcoxon(selected_diffs)
        print("HLA-{0}: W = {1}, total = {2}, P = {3} (wilcoxon)".format(locus, W, total, pval_wilcox))
        ax.set_xlim((-1.05*maxabsdiff, 1.05*maxabsdiff))
    ## additional formatting...
    fig.tight_layout()
    fig.savefig("compare.{0}.{1}.png".format(modelname1, modelname2), dpi=300, bbox_inches='tight')
