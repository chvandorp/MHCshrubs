"""@package mhcbind
An interface to NetMHCpan and other peptide-binding tools

The module contains fuctions for correctly calling NetMHCpan in parallel,
and reading the generated output from such a call.
The results can then be transformed into a MHC-similarity matrix.
"""
from __future__ import (print_function, division)
from builtins import str
import os
import multiprocessing
import subprocess
import numpy as np
import scipy.stats as sts
import csv
from mhcshrubs import auxiliary as aux
from mhcshrubs import (mhctools, progressbar)

supportedNetMHCpanVersions = ["2.8", "3.0"]
latestNetMHCpanVersion = "3.0" ## TODO


def mkNetMHCpanMhcStrings(mhcList, num=1, max_str_len=1024, sep=','):
    """
    Convert a list of MhcObjects into a string that can be passed to NetMHCpan.

    NetMHCpan only accepts strings that are shorter than 1024 characters.
    Therefore the string is split into shorter parts whenever the string would be too long.

    Args:
        mhcList (list): list of MhcObjects

    Kwargs:
        num (int): number of strings aimed for (default: 1)
        max_str_len (int): the maximal lengh of the string that may be passed
            to NetMHCpan (default: 1024)
        sep (str): the seperator for the alleles (default: ',')

    Returns:
        list of strings
    """
    def strLen(ml):
        return np.sum([len(mhc.NetMHC_str()) for mhc in ml]) + len(sep)*(len(ml) - 1)
    ## split the mhcList into num equal parts
    mhcLists = aux.split_list(mhcList, num)
    ## check that none of the resulting strings would be longer than max_str_len
    strLensOk = all([strLen(ml) < max_str_len for ml in mhcLists])
    if strLensOk:
        return [sep.join([mhc.NetMHC_str() for mhc in ml]) for ml in mhcLists]
    else: ## increase num and try again
        return mkNetMHCpanMhcStrings(mhcList, num=num+1, max_str_len=max_str_len, sep=sep)


def getNetMHCpanOutFileNames(infilename, cache_dir, mhc):
    rel_infilename = os.path.basename(infilename)
    namebase, ext = os.path.splitext(rel_infilename)
    ## TODO: check that extension is .fasta?
    outfilename = os.path.join(cache_dir, "{0}.{1}.xls".format(namebase, mhc.Newick_str()))
    errfilename = os.path.join(cache_dir, "{0}.{1}.err".format(namebase, mhc.Newick_str()))
    logfilename = os.path.join(cache_dir, "{0}.{1}.log".format(namebase, mhc.Newick_str()))
    return (outfilename, errfilename, logfilename)


def callNetMHCpan(jd):
    """
    Function for use with e.g. multiprocessing.Pool.map that correctly calls netmhcpan.

    Args:
        jd (dict): a dictionary containing
            - infilename: fasta file with proteins (or a list of peptides),
                must have extension '.fasta"
            - cache_dir: the directory to write the output to
            - mhc: MHC allele (MhcObject)
            - path: the path where NetMHCpan is located
            - version: the desired NetMHCpan version. Must be in supportedNetMHCpanVersions
            - dry_run: don't actually run netMHCpan
            - verbose: print messages

    Returns:
        A tuple consisting of
            - retcode: the result from subprocess.call
            - outfilename: the name of the file containing NetMHCpan results
    """
    assert(jd["version"] in supportedNetMHCpanVersions)
    retcode = 0 ## to be returned, also for dry-runs
    if jd["path"] is None:
        netMHCpanCmd = "netMHCpan"
    else: ## FIXME: this will only work at TBB
        netMHCpanCmd = os.path.join(jd["path"], "netMHCpan-{}".format(jd["version"]))
    outfilename, errfilename, logfilename = \
        getNetMHCpanOutFileNames(jd["infilename"], jd["cache_dir"], jd["mhc"])
    cmdlist = [netMHCpanCmd, ## netmchpan command
               "-f", jd["infilename"],
               "-xls" ,"-xlsfile", outfilename,
               "-a", jd["mhc"].NetMHC_str(),
               "-l", "9"] ## restrict to 9-mers (use "-p" for lists of peptides)
    if jd["verbose"]:
        message = "calling {0} with input {1}".format(cmdlist[0], cmdlist[2])
        if jd["dry_run"]: message = "[dry run] " + message
        retcode = subprocess.call(["echo", message])
    logfile = open(logfilename, "w")
    errfile = open(errfilename, "w")
    if not jd["dry_run"]:
        retcode = subprocess.call(cmdlist, stdout=logfile, stderr=errfile)
    ## close the logfile
    logfile.close()
    errfile.close()
    if jd["verbose"]:
        subprocess.call(["echo", "complete!"])
    return (retcode, outfilename)

def callNetMHCpanParallel(infilename, mhcList, path, cache_dir="/tmp/",
                          threads=1, version=latestNetMHCpanVersion,
                          verbose=False, dry_run=False):
    """
    split a list of MHC alleles in smaller chunks, and run callNetMHCpan in parallel,
    each with one chunk. return a list of results: the return codes and the outfilenames
    """
    ## make sure that the cache directory exists
    try:
        os.makedirs(cache_dir)
    except OSError:
        if not os.path.isdir(cache_dir): raise

    parallel = threads > 1

    jds = [{"infilename" : infilename,
            "cache_dir" : cache_dir,
            "mhc" : mhc,
            "path" : path,
            "version": version,
            "dry_run" : dry_run,
            "verbose" : False} for mhc in mhcList]
    if verbose:
        results = progressbar.parallel_map(callNetMHCpan, jds, max_workers=threads)
    else:
        with multiprocessing.Pool(processes=threads) as pool:
            results = pool.map(callNetMHCpan, jds)
    return results



def importNetMhcOutput(fileName, version="3.0"):
    """
    Import binding predictions produced by NetMHCpan and
    store them in a dictionary.
    """
    assert(version in supportedNetMHCpanVersions)
    fileHandle = open(fileName, 'r')
    csvReader = csv.reader(fileHandle, delimiter='\t', quotechar='\n')
    predictionTable = [line for line in csvReader]
    fileHandle.close()

    header1PredictionTable = predictionTable[0]
    header2PredictionTable = predictionTable[1]
    predictionTable = predictionTable[2:]

    ## TODO: infer form header
    if version == "2.8":
        fstMhcCol = 3
        numMhcCols = 3 ## aff, nM, rank
        peptideCol = 1
        posCol = 0
        proteinCol = 2
        affCol = 0
        nMCol = 1
        rankCol = 2
    elif version == "3.0":
        fstMhcCol = 3
        numMhcCols = 4 ## core, aff, nM, rank
        peptideCol = 1
        posCol = 0
        proteinCol = 2
        coreCol = 0
        affCol = 1
        nMCol = 2
        rankCol = 3

    mhc_names = [x for x in header1PredictionTable if x!='']
    predPerMhcDict = {}
    for row in predictionTable:
        for i, mhc_name in enumerate(mhc_names):
            ## put hla in proper format
            mhc = mhctools.MhcObject(mhctools.fromNetMhcFormat(mhc_name))
            predDict = {}
            pos = predDict["pos"] = int(row[posCol])
            pep = predDict["pep"] = row[peptideCol]
            prot = predDict["prot"] = row[proteinCol]
            pepid = "{0:s}_{1:d}_{2:s}".format(prot, pos, pep)
            predDict["aff"] = float(row[fstMhcCol + i*numMhcCols + affCol])
            predDict["nM"] = float(row[fstMhcCol + i*numMhcCols + nMCol])
            predDict["rank"] = float(row[fstMhcCol + i*numMhcCols + rankCol])
            if version == "3.0":
                predDict["core"] = row[fstMhcCol + i*numMhcCols + coreCol]
            if mhc in list(predPerMhcDict.keys()):
                predPerMhcDict[mhc][pepid] = predDict
            else:
                predPerMhcDict[mhc] = {pepid : predDict}
    return predPerMhcDict


def thresholdTransform(th):
    fun = lambda x: 1 if x <= th else 0
    return fun

def oneMinusLog50k(x):
    return 1.0 - np.log(x) / np.log(5e4)

def mkBindingCorrelationMatrix(predPerMhcDict, var="aff", fun=lambda x: x):
    """
    make a correlation mhc-correlation matrix based on binding affinity or rank.
    input:
    - predPerMhcDict -- binding predictions sorted by MHC and peptide id
    - var -- the key of the variable used for the correlations. Must be "aff", "nM", or "rank"
    - fun -- an optional transformation of the variable. identity by default.
      e.g. thresholdTransform(1.0) gives a 1% threshold cutoff
    output:
    - corrMat -- the desired correlation matrix
    - mhcList -- a list of MHCs in the right order
    """
    def calcMhcBindingCorr(pdd1, pdd2):
        meetpepids = sorted(list(set(pdd1.keys()).intersection(set(pdd2.keys()))))
        x1 = [fun(pdd1[pepid][var]) for pepid in meetpepids]
        x2 = [fun(pdd2[pepid][var]) for pepid in meetpepids]
        r, p = sts.pearsonr(x1, x2)
        return r

    mhcList = sorted(predPerMhcDict.keys())
    n = len(mhcList)
    corrMat = np.zeros((n, n))
    for i1, mhc1 in enumerate(mhcList):
        for i2, mhc2 in enumerate(mhcList):
            if i2 > i1:
                corrMat[i2, i1] = corrMat[i1, i2] = calcMhcBindingCorr(predPerMhcDict[mhc1], predPerMhcDict[mhc2])
            elif i1 == i2:
                corrMat[i1, i2] = 1.0

    return (corrMat, mhcList)
