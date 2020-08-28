"""
Defines a number of python versions of the JAGS models,
used for cross validation
"""
from __future__ import (print_function, division)
import scipy.stats as sts
import numpy as np
from mhcshrubs import papyjags as ppj
from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn


def normal_model(virus_load, censoring_code, virus_load_bound,
                 subject_idx, chain, hlaAlleles, allele_condition):
    """
    Compute the log posterior density (lpd) of an observation.
    The lpd can be conditioned on the patient expressing a particular allele.
    For instance when the allele data is (partially) missing.

    Args:
        virus_load: the VL of the patient.
        censoring_code: code indicating the kind of censoring on VL
        virus_load_bound: if censored, this is the boundary
        subject_idx: the index of the subject in the chain
        hlaAlleles: lists of HLA alleles, used only when allele_cond is not None
        allele_condition: if None, then every sample is used. If an allele is passed,
            then only samples where the subject expresses the allele are used.

    Returns:
        a tuple consisting of
            - the (conditional) lpd
            - the posterior probability of the condition
    """
    ## sample over chain
    states = ppj.transposeChain(chain)
    log_likes = []

    for i, state in enumerate(states):
        if allele_condition is not None:
            locus = allele_condition.locus
            keys = ["hla{0}{1}".format(locus, z) for z in [1,2]]
            ## NB: correct for indexing shift between JAGS and Python
            allele_idxs = [state[key][subject_idx]-1 for key in keys]
            alleles = [hlaAlleles[locus][idx] for idx in allele_idxs]
            ## continue for loop if allele is not expressed...
            if all([allele != allele_condition for allele in alleles]):
                continue
        tau_V = state["tau_V"] ## precision
        sigma_V = 1.0/np.sqrt(tau_V) ## precision is inverse variance
        alpha = state["alpha"]
        W = state["W"][subject_idx]
        ## compute log likelihood
        if censoring_code == defn.uncensored_code:
            log_likes.append(sts.norm.logpdf(virus_load, loc=W, scale=sigma_V))
        elif censoring_code == defn.left_censored_code: ## unknown virus load BELOW known value
            log_likes.append(sts.norm.logcdf(virus_load_bound, loc=W, scale=sigma_V))
        elif censoring_code == defn.right_censored_code: ## unknown virus load ABOVE known value
            ## scipy norm.logsf is the log-survival function (i.e. log(1-cdf))
            log_likes.append(sts.norm.logsf(virus_load_bound, loc=W, scale=sigma_V))
        elif censoring_code == aux.missing_code:
            log_likes.append(0.0)
    ## lpd = log(1/S sum p(y | theta)) with S the chain length
    lpd = aux.log_sum_exp(log_likes) - np.log(len(log_likes)) if len(log_likes) > 0 else np.nan
    pp_cond = len(log_likes) / len(states) if len(states) > 0 else np.nan
    return (lpd, pp_cond)
