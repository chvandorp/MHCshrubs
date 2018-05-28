"""@package crossval
Module for cross-validating the models.
"""
from __future__ import (print_function, division)
from builtins import zip
import multiprocessing
import time
import itertools
import numpy as np
import matplotlib.pyplot as plt
import csv
from mhcshrubs import (fittrees, mhcclus, mhctools, pymodels, sequences, progressbar, fetcher)
from mhcshrubs import auxiliary as aux
from mhcshrubs import papyjags as ppj

def crossValFun(task):
    """
    Auxiliary function for use with multiprocessing.Pool

    Args:
        task: is a dictionary containing the fields:
            - allele: the allele to be left out for the cross validation
            - tree: a newick tree giving the similarity structure of the HLA alleles
            - represDict: a dictionary of HLA representatives used in the tree
            - modelName: the name of the model
            - patientFileName: the name of the file with patient data
            - hlaFileName: the name of the file with allele frequency data
            - dry_run: don't sample the MCMC
            - use_cache: don't sample the MCMC, but use cached output
            - verbose: print progress
            - chain_len: the length of the MCMC
            - chain_thin: applied post sampling chain thinning
            - prior: string denoting the prior used for allele/node weights
            - multiplicity: string denoting the allele homozygosity model
            - hetr_adv: string (or None) denoting the heterozygote advantage

    Returns:
        dict: The result of fittrees.fitTreeWeights
    """
    if task["verbose"]:
        print("loo cross validation for", task["allele"], "with model", task["modelName"])
    result = fittrees.fitTreeWeights(task["dataDict"], task["parDict"], task["hlaFileName"], task["tree"],
                                     represDict=task["represDict"], chain_len=task["chain_len"],
                                     chain_thin=task["chain_thin"], leaveOutAlleles=[task["allele"]],
                                     modelName=task["modelName"], verbose=task["verbose"],
                                     render_figures=task["render_figures"],
                                     dry_run=task["dry_run"], use_cache=task["use_cache"],
                                     prior=task["prior"],
                                     multiplicity=task["multiplicity"], hetr_adv=task["hetr_adv"])
    return result


def posteriorProbAlleles(subject_idx, chain, hlaAlleles):
    """
    For a given subject, get the posterior probabilities of all alleles
    """
    pps = {locus : {allele : [0.0, 0.0] for allele in alleles}
           for locus, alleles in hlaAlleles.items()}
    for locus, alleles in hlaAlleles.items():
        for i in range(2):
            param_key = "hla{0}{1}".format(locus, i+1)
            allele_idxs = [x[subject_idx] for x in chain[param_key]]
            ## make sure to correct for the fact that JAGS starts indexing at 1
            counts = [sum(1 for idx in allele_idxs if j+1 == idx) for j, _ in enumerate(alleles)]
            probs = np.array(counts, dtype=float)/len(chain[param_key])
            for j, allele in enumerate(alleles):
                pps[locus][allele][i] = probs[j]
    return pps


def crossValidateAlleles(dataDict, parDict, hlaFileName, hlaNewickString, summaryFileName,
                         represDict=None, chain_len=1000, chain_thin=10,
                         modelName="hlatree", prior="norm", multiplicity="double", hetr_adv=None,
                         dry_run=False, use_cache=False, parallel=False, verbose=False):
    """
    Fit a model with an allele left out, then test the performance on the left out allele.

    Args:
        dataDict (dict): a dictionary with the following fields
            - TraitValues -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - CensBounds -- if left or right censored, this contains the bound
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO
        hlaFileName (str): name of file with HLA frequencies
        hlaNewickString (str): the Newick string representing the tree
        summaryFileName (str): name of file to write a summary of the cross validation to

    Kwargs:
        represDict (dict): a dictionary with representatives of alleles in the tree
        chain_len (int): the length of the MCMC (default: 1000)
        chain_thin (int): thinning applied to the MCMC to reduce autocorrelation
            (default: 10)
        prior (str): denoting prior for allele/node weights in the Bayesian model
            (default: "norm")
        multiplicity (str): string denoting the way allele homozygosity is implemented
            (default: "double")
        hetr_adv (str): string (or None) indicating the heterozygosity measure
            (default: None)
        dry_run (bool): don't run the sampler (default: False)
        parallel (bool): use multiple CPU cores (default: False)
        verbose (bool): print information about progress. This is rather limited for the
            parallel case (default: False)

    Returns:
        a dictionary mapping alleles to cross-validation results
    """
    ## load patient data (needed for validation below)
    patientLogVls = dataDict["TraitValues"] ## TODO rename
    VlCensCodes = dataDict["CensCodes"]
    VlBounds = dataDict["CensBounds"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())
    ## make a list of tasks
    tasks = [{"allele" : allele,
              "modelName" : modelName,
              "tree" : hlaNewickString,
              "represDict" : represDict,
              "dataDict" : dataDict,
              "parDict" : parDict,
              "hlaFileName" : hlaFileName,
              "chain_len" : chain_len,
              "chain_thin" : chain_thin,
              "render_figures" : False,
              "dry_run" : dry_run,
              "use_cache" : use_cache,
              "verbose" : False,
              "prior" : prior,
              "multiplicity" : multiplicity,
              "hetr_adv" : hetr_adv}
             for locus in loci for allele in hlaAlleles[locus]]
    ## decide the number of CPUs to use
    max_workers = max(1, multiprocessing.cpu_count()-1) if parallel else 1
    if verbose:
        print("running cross-validation on {} cpu cores".format(max_workers))
    results = progressbar.parallel_map(crossValFun, tasks, max_workers=max_workers)
    ## compute log-posterior densities and write to file
    with open(summaryFileName, 'w') as summaryFileHandle:
        ## for each left-out HLA allele and each model, predict VLs for left-out patients
        header = "allele\tmodel\tlpd\tallele_count\tsubjects"
        summaryFileHandle.write(header + "\n")
        for task, result in zip(tasks, results):
            chain = result["chain"]
            loPatientIdxs = result["leaveOutPatientIdxs"]
            allele = task["allele"]
            ## make a dictionary of patient index -> lpd
            lpds = {}
            ppes = {} ## posterior probability expression
            for idx in loPatientIdxs:
                ## compute subject-specific lpd
                trait = patientLogVls[idx]
                cens =  VlCensCodes[idx]
                bound = VlBounds[idx]
                clpd, ppc = pymodels.normal_model(trait, cens, bound, idx, chain,
                                                  hlaAlleles, allele)
                lpds[idx] = clpd
                ppes[idx] = ppc
            ## weighted sum of lpds. NB: when ppes[idx] == 0, lpds[ids] is NaN
            lpd = np.sum([lpds[idx]*ppes[idx] for idx in loPatientIdxs if ppes[idx] != 0.0])
            allele_count = np.sum([ppes[idx] for idx in loPatientIdxs])
            ## store the lpd and friends in the result dictionary
            result["lpd"] = lpd
            result["allele_count"] = allele_count
            result["lpds"] = lpds
            result["ppes"] = ppes
            subjectString = ";".join(["{0},{1},{2}".format(idx, lpds[idx], ppes[idx]) for idx in loPatientIdxs])
            summaryString = "{0}\t{1}\t{2}\t{3}\t{4}".format(allele, task["modelName"],
                                                        lpd, allele_count, subjectString)
            summaryFileHandle.write(summaryString + "\n")
    return dict(zip([task["allele"] for task in tasks], results))
