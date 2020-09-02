"""@package fittrees
Functions for fitting weights of the brances of a tree to disease trait data.

Todo:
    - move tree rendering (in fitTreeWeights) to another dedicated graphics function
    - make more plots of the input/output
    - new module for data imports? (including importPatientData etc.)
"""

from __future__ import (print_function, division)
from builtins import (map, str, zip, range)
import ete3
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
import csv
import warnings
from mhcshrubs import papyjags as ppj ## for JagsModel
from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn
from mhcshrubs import (statistics, stantools, colortrees, mhcclus, mhctools, fetcher, plots)

def getJagsCensCode(c):
    if c == defn.left_censored_code: return 0 ## upper bound
    elif c == defn.right_censored_code: return 1 ## lower bound
    else: return np.nan


def getStanCensCode(c):
    if c == defn.uncensored_code: return 0
    elif c == defn.left_censored_code: return 1
    elif c == defn.right_censored_code: return 2
    elif c == defn.missing_code: return 3
    else: raise Exception("unknown censoring code")

## functions to make JAGS-friendly data structures

def mkNodeMatrix(tree, nodeNames, hlaAlleles, represDict=None):
    """
    Transform a tree into a dictionary of matrices encoding the paths
    from the root to the leafs.

    Notice that an HLA allele could have no representative (e.g. Null alleles)
    and in that case, the row contains only zeros. The weight associated
    with a non-represented allele is therefore always zero. In the case of a
    Null allele, that is what we want (in most cases).

    Args:
        tree: an ete3 tree object
        nodeNames: a list of names of nodes in tree that must be used
            for the columns of the output matrices
        hlaAlleles: a dictionary (indexed by locus) of hla allele lists
            that must be used for the rows of the matrices

    Kwargs:
        represDict: hlaAlleles that are not in the tree must be represented
            by another allele. This relation is given by represDict.

    Returns:
        a dictionary (indexed by locus) of logical matrices.
    """
    paths = colortrees.getPathsFromEteTree(tree)
    nodeMatrix = {X : np.zeros((len(hlaAlleles[X]), len(nodeNames)))
                  for X in list(hlaAlleles.keys())}
    if represDict is not None:
        equivClassDict = aux.invertRepresDict(represDict)
    for (mhc_repr, path) in paths:
        mhc_repr = mhctools.MhcObject(mhc_repr) ## FIXME: tree should contain MhcObject
        if represDict is not None:
            ## notice that the represDict can contain more alleles than listed in hlaAlleles
            ## also note that representatives can have a different locus (e.g. in the KIR case)
            mhcs = [x for x in equivClassDict[mhc_repr] if x in aux.flatten(hlaAlleles.values())]
        else:
            mhcs = [mhc_repr]
        for mhc in mhcs:
            ## create the right row for each allele represented by mhc
            i = hlaAlleles[mhc.locus].index(mhc)
            for node in path:
                j = int(node.name)
                nodeMatrix[mhc.locus][i, j] = 1
    return nodeMatrix

def mkLengthVector(tree, nodeNames, lenfun=lambda x: x):
    """
    input: an ete tree and a list of node names (indices).
    output: a list of distances to the parent nodes in the tree,
    in the order of the list nodeNames.
    """
    return np.array([lenfun((tree&str(i)).dist) for i in nodeNames])


def uniqueHlaIndex(alleleVec):
    """
    input: binary vector indicating admissible HLA alleles
    at locus X and homolog i.
    output: np.nan if there are multiple admissible alleles,
    and the index of the HLA alleles (counting from 1) if there
    is exactly 1 possible allele.
    """
    indices = [i+1 for i,b in enumerate(alleleVec) if b]
    return np.nan if len(indices) != 1 else float(indices[0])

def mkAdmissibilityMatrix(patientAlleleVecs, X, i):
    """
    Args:
        patient allele vectors (Boolean), the locus and the homolog

    Returns:
        - a binary matrix with for each patient a row with ones when an
            allele is admissible (1).
        - derived from (1), a list of indices of known alleles, or np.nans
            for unknown alleles.
    """
    cXi = [[1 if b else 0 for b in a[i]] for a in patientAlleleVecs[X]]
    hlaXi = [uniqueHlaIndex(av) for av in cXi]
    return (np.array(cXi), np.array(hlaXi))

def mkAdmissibilityLists(patientAlleleVecs):
    """
    input: patient allele vectors (Boolean)
    output: lists of admissible alleles in order of locus, ploidy, subject
    """
    keys = sorted(patientAlleleVecs.keys())
    NumSubjects = len(patientAlleleVecs[keys[0]])
    Ploidy = len(patientAlleleVecs[keys[0]][0])
    def get_allele_idxs(bs):
        ## make a list of indices (starting at 1) from a list of booleans
        return [i+1 for i, b in enumerate(bs) if b]
    ## for each subject, make lists of admissible alleles
    alleleLists = [[[get_allele_idxs(patientAlleleVecs[X][s][i]) for X in keys]
                    for i in range(Ploidy)] for s in range(NumSubjects)]
    ## alleles are indexed 1, 2, ...
    return alleleLists


## functions to retrieve data from the chain

def getExpectedHlaCounts(chain, hlaAlleles, leaveOutPatientIdxs=[]):
    """
    Count the expected number of patients per HLA allele from MCMC output

    Todo:
        "heterozygous" and "homozygous" counts

    Args:
        chain -- MCMC output
        hlaAlleles -- dict of locus -> list of HLA allele names

    Kwargs:
        leaveOutPatientIdxs -- the indices of patients that
            must be exclude from the count. optional: if not specified,
            all patients are included

    Returns:
        dict of locus -> list of expectations
    """
    states = ppj.transposeChain(chain)
    loci = sorted(hlaAlleles.keys())
    counts = dict((X, [0.0 for _ in hlaAlleles[X]]) for X in loci)
    ## FIXME: better allele names in the JAGS model
    patientIdxs = [i for i in range(len(states[0]["hlaA1"])) if i not in leaveOutPatientIdxs]
    for state in states:
        for patient_idx in patientIdxs:
            for X in loci:
                for i in "12":
                    key = "hla{0}{1}".format(X, i)
                    hlaIdxXi = int(state[key][patient_idx])-1 ## JAGS counts from 1 (not 0)
                    counts[X][hlaIdxXi] += 1.0/len(states)
    return counts


def getExpectedHlaCountsStan(chain, hlaAlleles, leaveOutPatientIdxs=[]):
    """
    Count the expected number of patients per HLA allele from MCMC output

    Todo:
        "heterozygous" and "homozygous" counts

    Args:
        chain -- MCMC output
        hlaAlleles -- dict of locus -> list of HLA allele names

    Kwargs:
        leaveOutPatientIdxs -- the indices of patients that must be
            excluded from the count. optional: if not specified,
            all patients are included

    Returns:
        dict of locus -> list of expectations
    """
    loci = sorted(hlaAlleles.keys())
    counts = {X : [0.0 for _ in hlaAlleles[X]] for X in loci}
    ## TODO
    warnings.warn("fittrees.getExpectedHlaCountsStan is not implemented")
    return counts


## functions for the actual fitting...

def fitNullModel(dataDict, parDict, chain_len=1000, chain_thin=10, num_chains=4,
                 parallel=True, leaveOutAlleles=[], modelName="null",
                 verbose=False, dry_run=False, use_cache=False):
    """
    Fit the null model to the data.

    The null model is just a N(mu, sigma^2) model for the virus load with no effect of HLA.

    Args:
        dataDict -- a dictionary with the following fields
            - TraitValues -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - CensBounds -- if left or right censored, this contains the bound
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO

    Kwargs:
        chain_len -- the length of the MCMC (default: 1000)
        chain_thin -- the amount of thinning for the MCMC (default: 10)
        num_chains -- the number of independent chains (default: 4)
        parallel -- run independent chains concurrently (default: True)
        leaveOutAlleles -- mask VL of patients with a particular HLA allele (default: [])
        modelName -- identifier used for filenames (default: null)
        verbose -- print messages (default: False)
        dry_run -- don't run JAGS (default: False)

    Returns:
        a dictionary with
            - WAIC -- the WAIC of the run
            - chain -- a dictionary with the traces of the model parameters
            - leaveOutPatientIdxs -- indices of patients that were left out for cross validation
    """
    ## get the data
    traitValues = dataDict["TraitValues"]
    traitCensCodes = dataDict["CensCodes"]
    traitCensBounds = dataDict["CensBounds"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())

    leaveOutAlleleVec = {locus : [hla in leaveOutAlleles for hla in hlaAlleles[locus]] for locus in loci}

    ## find patients that should be left out
    def hasForbiddenAllele(av, locus):
        return any(aux.flatten([[a1 and a2 for a1, a2 in zip(av[i], leaveOutAlleleVec[locus])]
                                for i in range(2)]))

    leaveOutPatientIdxs = sorted(aux.flatten([[i for i, av in enumerate(patientAlleleVecs[X])
                                              if hasForbiddenAllele(av, X)] for X in loci]))

    ## set the traitValue (and censoring) of the left-out patients to NaN (and the proper censoring code)
    traitValues = [X if i not in leaveOutPatientIdxs else np.nan
                     for i, X in enumerate(traitValues)]
    traitCensCodes = [code if i not in leaveOutPatientIdxs else defn.missing_code
                   for i, code in enumerate(traitCensCodes)]
    traitCensBounds = [bound if i not in leaveOutPatientIdxs else defn.auxiliaryLowerCensBound
                       for i, bound in enumerate(traitCensBounds)]

    ## prepare data for JAGS
    N = len(traitValues)

    intervalCensorValues = list(map(getJagsCensCode, traitCensCodes))

    dataDict = {
        "N" : N,
        "V" : traitValues,
        "VLbound" : traitCensBounds,
        "intervalCensorValue" : intervalCensorValues
    }

    parameters = ["tau_V", "log_like", "V", "W", "R2"]

    ## add a specifier to the model name to indicate left out alleles
    if len(leaveOutAlleles) > 0:
        modelName += "." + ".".join(hla.Newick_str() for hla in leaveOutAlleles)

    ## choose the right jags model
    file_name = os.path.join(defn.ROOT_DIR, "jags/null-model.bug")

    ## choose a working folder
    work_folder = os.getcwd()
    ## check folder for output
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")
    ## TODO: make sure this folder and a figures subfolder exists

    path = os.path.join(outputFolder, "jags-cache")

    ## make a JAGS model object
    jm = ppj.JagsModel(file_name=file_name, model_name=modelName, path=path)

    ## run the JAGS model
    jm.sampling(data=dataDict, pars=parameters, iter=chain_len,
                chains=num_chains, warmup=chain_len, thin=chain_thin,
                verbose=verbose, dry_run=dry_run, parallel=parallel)
    if len(jm.sams) > 0:
        chain = ppj.mergeChains(jm.sams)
        ## calculate statistics...
        WAIC = statistics.calcWAIC(chain["log_like"], verbose)
        if len(jm.sams) > 1:
            Rhat = dict((pn, statistics.calcGelmanRubinRhat(jm.sams, pn)) for pn in parameters)
        else:
            Rhat = {}
    else:
        if verbose:
            print("no sample generated or no cached sample found.")
        chain = {}
        WAIC = np.nan
        Rhat = {}
    ## return a dictionary
    rd = {
        "WAIC" : WAIC,
        "Rhat" : Rhat,
        "chain" : chain,
        "leaveOutPatientIdxs" : leaveOutPatientIdxs
    }
    return rd


def fitTreeWeights(dataDict, parDict, hlaFileName, hlaNewickString,
                   represDict=None, chain_len=1000, chain_thin=10,
                   num_chains=4, parallel=False,
                   leaveOutAlleles=[], modelName="hlatree",
                   render_figures=True, verbose=False, dry_run=False,
                   use_cache=False, prior="norm",
                   multiplicity="double", hetr_adv=None):
    """
    Fit weights on the branches on an allele tree from a disease trait.

    Args:
        dataDict -- a dictionary with the following fields
            - Traits -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - CensBounds -- if left or right censored, this contains the bound
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO
        hlaFileName (str): name of file cointaining HLA frequencies
        hlaNewickString (str): the HLA tree in newick format

    Kwargs:
        represDict (dict):
            a dictionary specifying equivalence classes of HLA alleles
            i.e. patients HLA alleles are mapped to alleles in the tree. If None,
            it is assumed that all classes are of size 1 (default: None)
        chain_len (int): the length of the MCMC (default: 1000)
        chain_thin (int): the amount of thinning for the MCMC (default: 10)
        num_chains (int): the number of independent MCMCs
        parallel (bool): run multiple chains in parallel (or not)
        leaveOutAlleles (list):
            list of alleles that should not be used for the
            fitting, any patient with the allele is discarded (i.e. has missing VL).
        modelName (str): a string specifying a name for the model
        render_figures (bool): render trees, trace plots, violin plots etc.
        verbose (bool): print messages (default: False)
        dry_run (bool): don't run JAGS (default: False)
        prior (str):
            prior used for the branch weights choises: "norm", "dexp" (default: "norm")
        multiplicity (str):
            how to count weights. choices:
            - double (default): homozygous weights have a double effect
            - single: homozygous weights have a single effect
        hetr_adv (str):
            how (if) to count heterozygote advantage
            - None (default): don't add a parameter for the heterozygote advantage
            - discrete: look at discreet allele identity
            - tree: use the tree to determine the level of heterozygosity

    Returns:
        a dictionary containing:
            - tree: an ete3 object with weighted nodes
            - chain: the posterior distribution of the parameters and sampled
                HLA alleles
            - leaveOutPatientIdxs: indices of patients that were left out of the fit.
                If leaveOutAlleles is not empty
    """
    ## determine the data type
    categorical = True if "Categories" in dataDict.keys() else False
    ## get the data
    traitValues = dataDict["TraitValues"]
    traitCensCodes = dataDict["CensCodes"]
    traitCategories = dataDict["Categories"] if categorical else []
    traitCensBounds = [] if categorical else dataDict["CensBounds"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())

    if hlaFileName is not None:
        hlaCountDict = fetcher.importAlleleFrequencyData(hlaFileName)
    else: ## empty dictionary
        hlaCountDict = {}

    ## get indices of forbidden alleles
    leaveOutAlleleVec = {locus : [hla in leaveOutAlleles for hla in hlaAlleles[locus]]
                         for locus in loci}

    ## find patients that should be left out
    def hasForbiddenAllele(av, locus):
        return any(aux.flatten([[a1 and a2 for a1, a2 in zip(av[i], leaveOutAlleleVec[locus])]
                                for i in range(2)]))

    leaveOutPatientIdxs = sorted(aux.flatten([[i for i, av in enumerate(patientAlleleVecs[locus])
                                              if hasForbiddenAllele(av, locus)] for locus in loci]))

    ## set the VL (and censoring) of the left-out patients to NaN (and the proper censoring code)
    traitValues = [X if i not in leaveOutPatientIdxs else np.nan
                     for i, X in enumerate(traitValues)]
    traitCensCodes = [code if i not in leaveOutPatientIdxs else defn.missing_code
                   for i, code in enumerate(traitCensCodes)]
    traitCensBounds = [bound if i not in leaveOutPatientIdxs else defn.auxiliaryLowerCensBound
                for i, bound in enumerate(traitCensBounds)]

    ## make an ETE tree
    colorfun = lambda hlaStr: defn.locusColorDict[mhctools.MhcObject(hlaStr).locus]
    tree, nodeNames = colortrees.mkEteHlaTree(hlaNewickString, colorfun=colorfun)

    if verbose:
        print("number of nodes in tree:", len(nodeNames))

    ## prepare data for JAGS
    N = len(traitValues)
    H = {locus : len(hlaAlleles[locus]) for locus in loci}

    ## FIXME generic version for any locus
    cA1, hlaA1 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 0)
    cA2, hlaA2 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 1)
    cB1, hlaB1 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 0)
    cB2, hlaB2 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 1)
    cC1, hlaC1 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 0)
    cC2, hlaC2 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 1)

    nodeMatrix = mkNodeMatrix(tree, nodeNames, hlaAlleles, represDict=represDict)
    ## TESTING
    #fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    #ax1.pcolor(nodeMatrix["A"])
    #ax2.pcolor(nodeMatrix["B"])
    #ax3.pcolor(nodeMatrix["C"])
    #fig.savefig("test.png")
    ## END TESTING
    ## FIXME specify lenfun to transform distances
    lengthsEdgeToParent = mkLengthVector(tree, nodeNames)
    K = len(nodeNames)
    m = {locus : np.array([hlaCountDict[hla] if hla in list(hlaCountDict.keys()) else 0
                           for hla in hlaAlleles[locus]]) for locus in loci}
    M = {locus : np.sum(m[locus]) for locus in loci}
    intervalCensorValues = list(map(getJagsCensCode, traitCensCodes))

    ## FIXME generic loci
    dataDict = {"NodeMatrixA" : nodeMatrix['A'],
                "NodeMatrixB" : nodeMatrix['B'],
                "NodeMatrixC" : nodeMatrix['C'],
                "lengthsEdgeToParent" : lengthsEdgeToParent,
                "hlaA1" : hlaA1, "hlaA2" : hlaA2,
                "hlaB1" : hlaB1, "hlaB2" : hlaB2,
                "hlaC1" : hlaC1, "hlaC2" : hlaC2,
                "cA1" : cA1, "cA2" : cA2,
                "cB1" : cB1, "cB2" : cB2,
                "cC1" : cC1, "cC2" : cC2,
                "mA" : m['A'], "mB" : m['B'], "mC" : m['C'],
                "MA" : M['A'], "MB" : M['B'], "MC" : M['C'],
                "N" : N, "K" : K,
                "HA" : H['A'], "HB" : H['B'], "HC" : H['C'],
                "V" : traitValues,
                "intervalCensorValue" : intervalCensorValues,
                "VLbound" : traitCensBounds}

    parameters = ["betaA", "betaB", "betaC", "betaNodes", "meanHlaEffect",
                  "hlaA1", "hlaA2", "hlaB1", "hlaB2", "hlaC1", "hlaC2",
                  "pA", "pB", "pC", "tau_V", "tau_beta",
                  "log_like", "V", "W", "alpha", "R2"]
    if hetr_adv is not None:
        parameters += ["heterozygosity", "eta"]


    ## add a specifier to the model name to indicate left out alleles
    if len(leaveOutAlleles) > 0:
        modelName += "." + ".".join(hla.Newick_str() for hla in leaveOutAlleles)

    ## choose the right jags model
    work_folder = os.getcwd()
    if multiplicity == "double":
        if prior == "norm":
            file_name = os.path.join(defn.ROOT_DIR, "jags", "hla-tree-model.bug")
        elif prior == "dexp":
            file_name = os.path.join(defn.ROOT_DIR, "jags", "dexp-hla-tree-model.bug")
        else:
            raise Exception("invalid prior distribution")
    elif multiplicity == "single":
        if hetr_adv is None:
            if prior == "norm":
                file_name = os.path.join(defn.ROOT_DIR, "jags", "genhom-hla-tree-model.bug")
            elif prior == "dexp":
                file_name = os.path.join(defn.ROOT_DIR, "jags", "dexp-genhom-hla-tree-model.bug")
            else:
                raise Exception("invalid prior distribution")
        elif hetr_adv == "tree":
            if prior == "norm":
                file_name = os.path.join(defn.ROOT_DIR, "jags", "dosage-hla-tree-model.bug")
            elif prior == "dexp":
                file_name = os.path.join(defn.ROOT_DIR, "jags", "dexp-dosage-hla-tree-model.bug")
            else:
                raise Exception("invalid prior distribution")
        elif hetr_adv == "discrete":
            raise Exception("discrete heterozygosity measure not implemented (TODO)")
        else:
            raise Exception("invalid heterozygosity measure")
    else:
        raise Exception(f"multiplicity '{multiplicity}' not implemented")

    ## check folder for output
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")
    ## TODO: make sure this folder and a figures subfolder exists

    ## choose a working folder
    path = os.path.join(outputFolder, "jags-cache")

    ## make a JAGS model object
    jm = ppj.JagsModel(file_name=file_name, model_name=modelName, path=path, use_cache=use_cache)

    ## run the JAGS model
    if not use_cache:
        jm.sampling(data=dataDict, pars=parameters, iter=chain_len,
                    warmup=chain_len, thin=chain_thin, verbose=verbose,
                    dry_run=dry_run, chains=num_chains, parallel=parallel)
    if len(jm.sams) > 0:
        chain = ppj.mergeChains(jm.sams)
        ## compute statistics
        if len(jm.sams) > 1:
            Rhat = {pn : statistics.calcGelmanRubinRhat(jm.sams, pn) for pn in parameters}
        else: ## Rhat is not well defined when there is only one chain
            Rhat = {}
        WAIC = statistics.calcWAIC(chain["log_like"], verbose)
        ## Bayesian R^2
        R2 = (np.mean(chain["R2"]), np.percentile(chain["R2"], 2.5), np.percentile(chain["R2"], 97.5))
        if verbose: print(f"R^2 = {R2[0]:0.3f} (95% CrI: [{R2[1]:0.3f},{R2[2]:0.3f}])")
        ## add the estimates to the tree
        betaNodes = [np.array(list(map(float, b))) for b in aux.transpose(chain["betaNodes"])]
        ## TODO: float conversion no longer needed here
        colortrees.addNodeFeatures(tree, nodeNames, "beta", betaNodes)
        colortrees.addCumulativeFeature(tree, "beta", scale_dist=True)
        expectedHlaCounts = getExpectedHlaCounts(chain, hlaAlleles, leaveOutPatientIdxs)
        colortrees.addLeafBars(tree, hlaAlleles, expectedHlaCounts, represDict=represDict)
        if hetr_adv is not None:
            etas = chain["eta"]
            print("heterozygous advantage eta. mean:", np.mean(etas), "sd:", np.std(etas))
            heterozygosities = aux.flatten(chain["heterozygosity"])
            m_hetr = np.mean(heterozygosities)
            l_hetr = np.percentile(heterozygosities, 2.5)
            u_hetr = np.percentile(heterozygosities, 97.5)
            print(f"heterozygosity mean: {m_hetr:0.2f}, 2.5-97.5 percentile: [{l_hetr:0.2f}, {u_hetr:0.2f}]")
    else:
        if verbose:
            warnings.warn("warning: no sample generated or no cached sample found.")
        chain = {}
        Rhat = {}
        WAIC = np.nan
        R2 = (np.nan, np.nan, np.nan)

    ## draw tree with colors indicating weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder, "figures",
        f"colored-mhc-tree.{modelName}.png")
    colortrees.colorEteTree(tree, cfname="beta", wfname="beta", cffun=np.mean,
                            wffun=aux.getMass, cf_scale_dist=True)
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder, "figures",
        f"colored-mhc-utree.{modelName}.png")
    if render_figures:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## draw tree with colors indicating cumulative weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder, "figures",
        f"colored-mhc-tree.{modelName}.cumulative.png")
    legend_fig, (cax, wax) = plt.subplots(2, 1, figsize=(3, 2))
    colortrees.colorEteTree(tree, cfname="beta_cumulative", wfname="beta",
                            cffun=np.mean, wffun=(lambda x: abs(np.mean(x))),
                            cf_scale_dist=False, wf_scale_dist=True, cax=cax, wax=wax)
    legend_fig.tight_layout()
    legend_filename = os.path.join(outputFolder, "figures" , f"legend.{modelName}.png")
    if render_figures:
        legend_fig.savefig(legend_filename, dpi=300, bbox_inches='tight')
    plt.close(legend_fig)
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder, "figures",
        f"colored-mhc-utree.{modelName}.cumulative.png")
    if render_figures:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
            ## make a couple of trace plots
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## make trace plots
    if render_figures:
        if aux.isXServerAvailable():
            figFileNameTP = os.path.join(outputFolder,
              "figures/trace-plots.{}.png".format(modelName))
            parametersTP = [
                ("alpha", None), ## TODO translate parameter names
                ("tau_V", None),
                ("tau_beta", None),
                ("betaNodes", np.random.randint(0, K)),
                ("betaNodes", np.random.randint(0, K)),
                ("betaNodes", np.random.randint(0, K)),
                ("betaA", np.random.randint(0, H["A"])),
                ("betaB", np.random.randint(0, H["B"])),
                ("betaC", np.random.randint(0, H["C"])),
                ("pA", np.random.randint(0, H["A"])),
                ("pB", np.random.randint(0, H["B"])),
                ("pC", np.random.randint(0, H["C"]))
            ]
            statistics.mkTracePlots(figFileNameTP, jm.sams, parametersTP, Rhat)
            ## allele weight plots
            figFileNameWeights = os.path.join(outputFolder,
              "figures/weight-plots.{}.png".format(modelName))
            plots.mkAlleleWeightPlots(figFileNameWeights, chain, hlaAlleles)
            ## allele weight plots
            figFileNameFreqs = os.path.join(outputFolder, "figures",
                f"freq-plots.{modelName}.png")
            plots.mkAlleleFreqPlots(figFileNameFreqs, chain, hlaAlleles, m)
        elif verbose:
            warnings.warn("no X-server avaliable. Trace plot not rendered.")
    ## return results
    rd = { ## result dictionary
        "WAIC" : WAIC,
        "Rhat" : Rhat,
        "R2" : R2,
        #"tree" : tree, ## FIXME: concurrency needs to pickle stuff: confilct with ete3.Tree?
        "chain" : chain,
        "leaveOutPatientIdxs" : leaveOutPatientIdxs
    }
    return rd


def fitTreeWeightsCat(dataDict, parDict, hlaFileName, hlaNewickString,
                   represDict=None, chain_len=1000, chain_thin=10,
                   num_chains=4, parallel=False,
                   leaveOutAlleles=[], modelName="hlatree",
                   render_figures=True, verbose=False, dry_run=False,
                   use_cache=False, prior="norm",
                   multiplicity="double", hetr_adv=None):
    """
    Fit weights on the branches on an allele tree from a disease trait.
    Specialized method for categorical data.
    @todo: make a single functio with options?

    Args:
        dataDict -- a dictionary with the following fields
            - Traits -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - Categories -- list of categories
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO
        hlaFileName (str): name of file cointaining HLA frequencies
        hlaNewickString (str): the HLA tree in newick format

    Kwargs:
        represDict (dict):
            a dictionary specifying equivalence classes of HLA alleles
            i.e. patients HLA alleles are mapped to alleles in the tree. If None,
            it is assumed that all classes are of size 1 (default: None)
        chain_len (int): the length of the MCMC (default: 1000)
        chain_thin (int): the amount of thinning for the MCMC (default: 10)
        num_chains (int): the number of independent MCMCs
        parallel (bool): run multiple chains in parallel (or not)
        leaveOutAlleles (list):
            list of alleles that should not be used for the
            fitting, any patient with the allele is discarded (i.e. has missing VL).
        modelName (str): a string specifying a name for the model
        render_figures (bool): render trees, trace plots, violin plots etc.
        verbose (bool): print messages (default: False)
        dry_run (bool): don't run JAGS (default: False)
        prior (str):
            prior used for the branch weights choises: "norm", "dexp" (default: "norm")
        multiplicity (str):
            how to count weights. choices:
            - double (default): homozygous weights have a double effect
            - single: homozygous weights have a single effect
        hetr_adv (str):
            how (if) to count heterozygote advantage
            - None (default): don't add a parameter for the heterozygote advantage
            - discrete: look at discreet allele identity
            - tree: use the tree to determine the level of heterozygosity

    Returns:
        a dictionary containing:
            - tree: an ete3 object with weighted nodes
            - chain: the posterior distribution of the parameters and sampled
                HLA alleles
            - leaveOutPatientIdxs: indices of patients that were left out of the fit.
                If leaveOutAlleles is not empty
    """
    ## get the data
    traitValues = dataDict["TraitValues"]
    traitCensCodes = dataDict["CensCodes"]
    traitCategories = dataDict["Categories"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())

    if hlaFileName is not None:
        hlaCountDict = fetcher.importAlleleFrequencyData(hlaFileName)
    else: ## empty dictionary
        hlaCountDict = {}

    ## get indices of forbidden alleles
    leaveOutAlleleVec = {locus : [hla in leaveOutAlleles for hla in hlaAlleles[locus]]
                         for locus in loci}

    ## find patients that should be left out
    def hasForbiddenAllele(av, locus):
        return any(aux.flatten([[a1 and a2 for a1, a2 in zip(av[i], leaveOutAlleleVec[locus])]
                                for i in range(2)]))

    leaveOutPatientIdxs = sorted(aux.flatten([[i for i, av in enumerate(patientAlleleVecs[locus])
                                              if hasForbiddenAllele(av, locus)] for locus in loci]))

    ## set the VL (and censoring) of the left-out patients to NaN (and the proper censoring code)
    traitValues = [X if i not in leaveOutPatientIdxs else np.nan
                     for i, X in enumerate(traitValues)]
    traitCensCodes = [code if i not in leaveOutPatientIdxs else defn.missing_code
                   for i, code in enumerate(traitCensCodes)]

    ## mape trait-values to 0/1
    if len(traitCategories) != 2:
        raise Exception("invalid number of categories: only binary outcomes are implemented")
    def cat_map(x):
        if x == traitCategories[0]:
            return 0
        elif x == traitCategories[1]:
            return 1
        raise Exception(f"trait value {x} not in category list")
    traitValues = [cat_map(x) if c == defn.uncensored_code else np.nan
                   for x, c in zip(traitValues, traitCensCodes)]

    ## make an ETE tree
    colorfun = lambda hlaStr: defn.locusColorDict[mhctools.MhcObject(hlaStr).locus]
    tree, nodeNames = colortrees.mkEteHlaTree(hlaNewickString, colorfun=colorfun)

    if verbose:
        print("number of nodes in tree:", len(nodeNames))

    ## prepare data for JAGS
    N = len(traitValues)
    H = {locus : len(hlaAlleles[locus]) for locus in loci}

    ## FIXME generic version for any locus
    cA1, hlaA1 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 0)
    cA2, hlaA2 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 1)
    cB1, hlaB1 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 0)
    cB2, hlaB2 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 1)
    cC1, hlaC1 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 0)
    cC2, hlaC2 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 1)

    nodeMatrix = mkNodeMatrix(tree, nodeNames, hlaAlleles, represDict=represDict)
    ## TESTING
    #fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    #ax1.pcolor(nodeMatrix["A"])
    #ax2.pcolor(nodeMatrix["B"])
    #ax3.pcolor(nodeMatrix["C"])
    #fig.savefig("test.png")
    ## END TESTING
    ## FIXME specify lenfun to transform distances
    lengthsEdgeToParent = mkLengthVector(tree, nodeNames)
    K = len(nodeNames)
    m = {locus : np.array([hlaCountDict[hla] if hla in list(hlaCountDict.keys()) else 0
                           for hla in hlaAlleles[locus]]) for locus in loci}
    M = {locus : np.sum(m[locus]) for locus in loci}

    ## FIXME generic loci
    dataDict = {"NodeMatrixA" : nodeMatrix['A'],
                "NodeMatrixB" : nodeMatrix['B'],
                "NodeMatrixC" : nodeMatrix['C'],
                "lengthsEdgeToParent" : lengthsEdgeToParent,
                "hlaA1" : hlaA1, "hlaA2" : hlaA2,
                "hlaB1" : hlaB1, "hlaB2" : hlaB2,
                "hlaC1" : hlaC1, "hlaC2" : hlaC2,
                "cA1" : cA1, "cA2" : cA2,
                "cB1" : cB1, "cB2" : cB2,
                "cC1" : cC1, "cC2" : cC2,
                "mA" : m['A'], "mB" : m['B'], "mC" : m['C'],
                "MA" : M['A'], "MB" : M['B'], "MC" : M['C'],
                "N" : N, "K" : K,
                "HA" : H['A'], "HB" : H['B'], "HC" : H['C'],
                "Y" : traitValues}

    parameters = ["betaA", "betaB", "betaC", "betaNodes", "meanHlaEffect",
                  "hlaA1", "hlaA2", "hlaB1", "hlaB2", "hlaC1", "hlaC2",
                  "pA", "pB", "pC", "tau_beta",
                  "log_like", "Y", "Q", "alpha", "R2"]
    if hetr_adv is not None:
        parameters += ["heterozygosity", "eta"]


    ## add a specifier to the model name to indicate left out alleles
    if len(leaveOutAlleles) > 0:
        modelName += "." + ".".join(hla.Newick_str() for hla in leaveOutAlleles)

    ## choose the right jags model
    work_folder = os.getcwd()
    if multiplicity == "double":
        if prior == "norm":
            file_name = os.path.join(defn.ROOT_DIR, "jags", "hla-tree-model-binary.bug")
        elif prior == "dexp":
            file_name = os.path.join(defn.ROOT_DIR, "jags", "dexp-hla-tree-model-binary.bug")
        else:
            raise Exception("invalid prior distribution")
    else:
        ## @todo: implement more models
        raise Exception(f"multiplicity '{multiplicity}' not implemented")

    ## check folder for output
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")
    ## TODO: make sure this folder and a figures subfolder exists

    ## choose a working folder
    path = os.path.join(outputFolder, "jags-cache")

    ## make a JAGS model object
    jm = ppj.JagsModel(file_name=file_name, model_name=modelName, path=path, use_cache=use_cache)

    ## run the JAGS model
    if not use_cache:
        jm.sampling(data=dataDict, pars=parameters, iter=chain_len,
                    warmup=chain_len, thin=chain_thin, verbose=verbose,
                    dry_run=dry_run, chains=num_chains, parallel=parallel)
    if len(jm.sams) > 0:
        chain = ppj.mergeChains(jm.sams)
        ## compute statistics
        if len(jm.sams) > 1:
            Rhat = {pn : statistics.calcGelmanRubinRhat(jm.sams, pn) for pn in parameters}
        else: ## Rhat is not well defined when there is only one chain
            Rhat = {}
        WAIC = statistics.calcWAIC(chain["log_like"], verbose)
        ## Bayesian R^2
        R2 = (np.mean(chain["R2"]), np.percentile(chain["R2"], 2.5), np.percentile(chain["R2"], 97.5))
        if verbose: print(f"R^2 = {R2[0]:0.3f} (95% CrI: [{R2[1]:0.3f},{R2[2]:0.3f}])")
        ## add the estimates to the tree
        betaNodes = [np.array(list(map(float, b))) for b in aux.transpose(chain["betaNodes"])]
        ## TODO: float conversion no longer needed here
        colortrees.addNodeFeatures(tree, nodeNames, "beta", betaNodes)
        colortrees.addCumulativeFeature(tree, "beta", scale_dist=True)
        expectedHlaCounts = getExpectedHlaCounts(chain, hlaAlleles, leaveOutPatientIdxs)
        colortrees.addLeafBars(tree, hlaAlleles, expectedHlaCounts, represDict=represDict)
        if hetr_adv is not None:
            etas = chain["eta"]
            print("heterozygous advantage eta. mean:", np.mean(etas), "sd:", np.std(etas))
            heterozygosities = aux.flatten(chain["heterozygosity"])
            m_hetr = np.mean(heterozygosities)
            l_hetr = np.percentile(heterozygosities, 2.5)
            u_hetr = np.percentile(heterozygosities, 97.5)
            print(f"heterozygosity mean: {m_hetr:0.2f}, 2.5-97.5 percentile: [{l_hetr:0.2f}, {u_hetr:0.2f}]")
    else:
        if verbose:
            warnings.warn("warning: no sample generated or no cached sample found.")
        chain = {}
        Rhat = {}
        WAIC = np.nan
        R2 = (np.nan, np.nan, np.nan)

    ## draw tree with colors indicating weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder, "figures",
        f"colored-mhc-tree.{modelName}.png")
    colortrees.colorEteTree(tree, cfname="beta", wfname="beta", cffun=np.mean,
                            wffun=aux.getMass, cf_scale_dist=True)
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder, "figures",
        f"colored-mhc-utree.{modelName}.png")
    if render_figures:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## draw tree with colors indicating cumulative weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder, "figures",
        f"colored-mhc-tree.{modelName}.cumulative.png")
    legend_fig, (cax, wax) = plt.subplots(2, 1, figsize=(3, 2))
    colortrees.colorEteTree(tree, cfname="beta_cumulative", wfname="beta",
                            cffun=np.mean, wffun=(lambda x: abs(np.mean(x))),
                            cf_scale_dist=False, wf_scale_dist=True, cax=cax, wax=wax)
    legend_fig.tight_layout()
    legend_filename = os.path.join(outputFolder, "figures" , f"legend.{modelName}.png")
    if render_figures:
        legend_fig.savefig(legend_filename, dpi=300, bbox_inches='tight')
    plt.close(legend_fig)
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder, "figures",
        f"colored-mhc-utree.{modelName}.cumulative.png")
    if render_figures:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
            ## make a couple of trace plots
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## make trace plots
    if render_figures:
        if aux.isXServerAvailable():
            figFileNameTP = os.path.join(outputFolder,
              "figures/trace-plots.{}.png".format(modelName))
            parametersTP = [
                ("alpha", None), ## TODO translate parameter names
                ("tau_beta", None),
                ("betaNodes", np.random.randint(0, K)),
                ("betaNodes", np.random.randint(0, K)),
                ("betaNodes", np.random.randint(0, K)),
                ("betaA", np.random.randint(0, H["A"])),
                ("betaB", np.random.randint(0, H["B"])),
                ("betaC", np.random.randint(0, H["C"])),
                ("pA", np.random.randint(0, H["A"])),
                ("pB", np.random.randint(0, H["B"])),
                ("pC", np.random.randint(0, H["C"]))
            ]
            statistics.mkTracePlots(figFileNameTP, jm.sams, parametersTP, Rhat)
            ## allele weight plots
            figFileNameWeights = os.path.join(outputFolder,
              "figures/weight-plots.{}.png".format(modelName))
            plots.mkAlleleWeightPlots(figFileNameWeights, chain, hlaAlleles)
            ## allele weight plots
            figFileNameFreqs = os.path.join(outputFolder, "figures",
                f"freq-plots.{modelName}.png")
            plots.mkAlleleFreqPlots(figFileNameFreqs, chain, hlaAlleles, m)
        elif verbose:
            warnings.warn("no X-server avaliable. Trace plot not rendered.")
    ## return results
    rd = { ## result dictionary
        "WAIC" : WAIC,
        "Rhat" : Rhat,
        "R2" : R2,
        #"tree" : tree, ## FIXME: concurrency needs to pickle stuff: confilct with ete3.Tree?
        "chain" : chain,
        "leaveOutPatientIdxs" : leaveOutPatientIdxs
    }
    return rd




def fitPMMweights(dataDict, parDict, hlaFileName, hlaCovMat, hlaCovMatHeader,
                  represDict=None, chain_len=1000, chain_thin=10, num_chains=4,
                  parallel=True, leaveOutAlleles=[], modelName="pmm",
                  render_figures=True, verbose=False, dry_run=False, use_cache=False):
    """
    Use the PMM (phylogenetic mixed model) to find HLA associations.

    A distance matrix defining dissimilarities between HLA alleles is
    used to define a correlation structure between the alleles.
    Then, a JAGS model is used to find weights.

    Args:
        dataDict -- a dictionary with the following fields
            - Traits -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - CensBounds -- if left or right censored, this contains the bound
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO
        hlaFileName -- name of file containing HLA frequencies
        hlaCovMat -- similarities between HLA alleles (should be positive definite)
        hlaCovMatHeader -- order of alleles used in the hlaCovMat

    Kwargs:
        represDict -- representatives of alleles with that prefectly covary
        chain_len -- the length of the MCMC (default: 1000)
        chain_thin -- the amount of thinning for the MCMC (default: 10)
        num_chains -- number of independent runs (default: 4)
        parallel -- run independent runs concurrently (default: True)
        leaveOutAlleles -- list of alleles that should not be used for the
            fitting, any patient with the allele is discarded (i.e. has missing VL).
        modelName -- a string specifying a name for the model
        render_figures -- if True, render trace plots, violin plots, etc.
        verbose -- print messages (default: False)
        dry_run -- don't run JAGS (default: False)

    Returns:
        a dictionary containing:
            - WAIC -- the WAIC of the model
            - chain -- the posterior distribution of the parameters and sampled
                HLA alleles
            - leaveOutPatientIdxs -- indices of patients that were left out of
                the fit. if leaveOutAlleles is not empty
    """
    ## get the data
    traitValues = dataDict["TraitValues"]
    traitCensCodes = dataDict["CensCodes"]
    traitCensBounds = dataDict["CensBounds"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())

    if hlaFileName is not None:
        hlaCountDict = fetcher.importAlleleFrequencyData(hlaFileName)
    else: ## empty dictionary
        hlaCountDict = {}

    ## get indices of forbidden alleles
    leaveOutAlleleVec = {locus : [hla in leaveOutAlleles for hla in hlaAlleles[locus]] for locus in loci}

    ## find patients that should be left out
    def hasForbiddenAllele(av, locus):
        return any(aux.flatten([[a1 and a2 for a1, a2 in zip(av[i], leaveOutAlleleVec[locus])]
                                for i in range(2)]))

    leaveOutPatientIdxs = sorted(aux.flatten([[i for i, av in enumerate(patientAlleleVecs[locus])
                                              if hasForbiddenAllele(av, locus)] for locus in loci]))

    ## set the VL (and censoring) of the left-out patients to NaN (and the proper censoring code)
    traitValues = [vl if i not in leaveOutPatientIdxs else np.nan
                     for i, vl in enumerate(traitValues)]
    traitCensCodes = [code if i not in leaveOutPatientIdxs else defn.missing_code
                   for i, code in enumerate(traitCensCodes)]
    traitCensBounds = [bound if i not in leaveOutPatientIdxs else defn.auxiliaryLowerCensBound
                for i, bound in enumerate(traitCensBounds)]

    ## test that the covariance matrix is positive definite
    if not aux.is_pos_def(hlaCovMat):
        raise Exception("Non-positive definite covariance matrix")

    if verbose:
        print("number of HLA representatives: {}".format(hlaCovMat.shape[0]))
    ## prepare data for JAGS
    N = len(traitValues)
    H = {locus : len(hlaAlleles[locus]) for locus in loci}

    cA1, hlaA1 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 0)
    cA2, hlaA2 = mkAdmissibilityMatrix(patientAlleleVecs, "A", 1)
    cB1, hlaB1 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 0)
    cB2, hlaB2 = mkAdmissibilityMatrix(patientAlleleVecs, "B", 1)
    cC1, hlaC1 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 0)
    cC2, hlaC2 = mkAdmissibilityMatrix(patientAlleleVecs, "C", 1)

    m = {locus : np.array([hlaCountDict[hla]
                 if hla in list(hlaCountDict.keys()) else 0
                 for hla in hlaAlleles[locus]])
             for locus in loci}
    M = {locus : np.sum(m[locus]) for locus in loci}

    intervalCensorValues = list(map(getJagsCensCode, traitCensCodes))

    ## make a translation from HLA alleles to position in hlaCovMat (Sigma)
    ## notice the "+1" to correct for JAGS indexing of arrays
    represList = {locus : [hlaCovMatHeader.index(represDict[allele]) + 1
                           for allele in hlaAlleles[locus]]
                  for locus in loci}

    dataDict = {"Sigma" : hlaCovMat,
                "R" : len(hlaCovMatHeader),
                "IdxsSigmaA" : represList['A'],
                "IdxsSigmaB" : represList['B'],
                "IdxsSigmaC" : represList['C'],
                "hlaA1" : hlaA1, "hlaA2" : hlaA2,
                "hlaB1" : hlaB1, "hlaB2" : hlaB2,
                "hlaC1" : hlaC1, "hlaC2" : hlaC2,
                "cA1" : cA1, "cA2" : cA2,
                "cB1" : cB1, "cB2" : cB2,
                "cC1" : cC1, "cC2" : cC2,
                "mA" : m['A'], "mB" : m['B'], "mC" : m['C'],
                "MA" : M['A'], "MB" : M['B'], "MC" : M['C'],
                "HA" : H['A'], "HB" : H['B'], "HC" : H['C'],
                "N" : N, "V" : traitValues,
                "intervalCensorValue" : intervalCensorValues,
                "VLbound" : traitCensBounds}

    parameters = ["alpha", "betaA", "betaB", "betaC", "meanHlaEffect",
                  "hlaA1", "hlaA2", "hlaB1", "hlaB2", "hlaC1", "hlaC2",
                  "pA", "pB", "pC", "tau_V", "tau_beta",
                  "log_like", "V", "W", "R2"]

    ## add a specifier to the model name to indicate left out alleles
    if len(leaveOutAlleles) > 0:
        modelName += "." + ".".join(hla.Newick_str() for hla in leaveOutAlleles)

    work_folder = os.getcwd()
    ## check folder for output
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")

    ## choose the right jags model
    file_name = os.path.join(defn.ROOT_DIR, "jags", "hla-covar-model.bug")

    ## choose a working folder
    path = os.path.join(outputFolder, "jags-cache")

    ## make a JAGS model object
    jm = ppj.JagsModel(file_name=file_name, model_name=modelName,
                       path=path, use_cache=use_cache)

    ## run the JAGS model
    if not use_cache:
        jm.sampling(data=dataDict, pars=parameters, iter=chain_len,
                    warmup=chain_len, thin=chain_thin, chains=num_chains,
                    verbose=verbose, dry_run=dry_run, parallel=parallel)

    if len(jm.sams) > 0:
        chain = ppj.mergeChains(jm.sams)
        WAIC = statistics.calcWAIC(chain["log_like"], verbose=verbose)
        if len(jm.sams) > 1:
            Rhat = {pn : statistics.calcGelmanRubinRhat(jm.sams, pn) for pn in parameters}
        else: ## Rhat is not well defined when there is only one chain
            Rhat = {}
        ## Bayesian R^2
        R2 = (np.mean(chain["R2"]), np.percentile(chain["R2"], 2.5), np.percentile(chain["R2"], 97.5))
        if verbose: print(f"R^2 = {R2[0]:0.3f} (95% CrI: [{R2[1]:0.3f},{R2[2]:0.3f}])")
    else:
        if verbose:
            print("no sample generated or no cached sample found.")
        chain = {}
        WAIC = np.nan
        Rhat = {}
        R2 = (np.nan, np.nan, np.nan)

    if render_figures:
        if aux.isXServerAvailable():
            figFileNameTP = os.path.join(outputFolder,
              "figures", "trace-plots.{}.png".format(modelName))
            parametersTP = [
                ("alpha", None), ## TODO: translate parameter names
                ("tau_V", None),
                ("tau_beta", None),
                ("betaA", np.random.randint(0, H["A"])),
                ("betaB", np.random.randint(0, H["B"])),
                ("betaC", np.random.randint(0, H["C"])),
                ("pA", np.random.randint(0, H["A"])),
                ("pB", np.random.randint(0, H["B"])),
                ("pC", np.random.randint(0, H["C"]))
            ]
            statistics.mkTracePlots(figFileNameTP, jm.sams, parametersTP, Rhat)
            ## allele weight plots
            figFileNameWeights = os.path.join(outputFolder,
              "figures", f"weight-plots.{modelName}.png")
            plots.mkAlleleWeightPlots(figFileNameWeights, chain, hlaAlleles)
            ## allele weight plots
            figFileNameFreqs = os.path.join(outputFolder,
              "figures", f"freq-plots.{modelName}.png")
            plots.mkAlleleFreqPlots(figFileNameFreqs, chain, hlaAlleles, m)
        elif verbose:
            warnings.warn("No X-server avaliable. Figures not rendered.")

    ## return the results
    rd = {
        "chain" : chain,
        "WAIC" : WAIC,
        "R2" : R2,
        "Rhat" : Rhat,
        "leaveOutPatientIdxs" : leaveOutPatientIdxs
    }


def fitTreeWeightsStan(dataDict, parDict, hlaFileName, hlaNewickString,
                       represDict=None, chain_len=1000, chain_thin=1, num_chains=4,
                       leaveOutAlleles=[], modelName="hlatree",
                       render_trees=True, verbose=False, dry_run=False, use_cache=False,
                       prior="norm", multiplicity="double", hetr_adv=None, wbic_sampling=True):
    """
    Fit the weights of the branches of a tree to allele-trait data using Stan.

    Args:
        dataDict -- a dictionary with the following fields
            - TraitValues -- the response variable
            - CensCodes -- codes indicating censoring of the response variable
            - CensBounds -- if left or right censored, this contains the bound
            - AlleleVecs -- vectors indicating alleles
            - Alleles -- a list of alleles in the correct order
        parDict -- a dictionary with the following fields
            - outputFolder -- folder for writing output
            - TODO
        hlaFileName -- name of file cointaining HLA frequencies
        hlaNewickString -- the HLA tree in newick format

    Kwargs:
        represDcit -- a dictionary specifying equivalence classes of HLA alleles
            i.e. patients HLA alleles are mapped to alleles in the tree. In None,
            it is assumed that all classes are of size 1 (default: None)
        chain_len -- the length of the MCMC (default: 1000)
        chain_thin -- the amount of thinning for the MCMC (default: 10)
        num_chains -- number of independent runs
        leaveOutAlleles -- list of alleles that should not be used for the
            fitting, any patient with the allele is discarded (i.e. has missing VL).
        modelName -- a string specifying a name for the model
        verbose -- print messages (default: False)
        dry_run -- don't run Stan (default: False)
        prior -- prior used for the branch weights choises: "norm", "dexp"
            (default: "norm")
        multiplicity -- how to count weights. choices:
            - double: homozygous weights have a double effect
            - single: homozygous weights have a single effect
            (default: double)
        hetr_adv (str):
            how (if) to count heterozygote advantage
            - None (default): don't add a parameter for the heterozygote advantage
            - discreet: look at discreet allele identity
            - tree: use the tree to determine the level of heterozygosity
        wbic_sampling -- if True (default) also sample at a higher temperature to
            calculate WBIC

    Returns:
        a dictionary containing:
            - tree -- an ete3 object with weighted nodes
            - chain -- the posterior distribution of the parameters and sampled
                HLA alleles
            - leaveOutPatientIdxs -- indices of patients that were left out of
                the fit. if leaveOutAlleles is not empty
    """
    ## get the data
    traitValues = dataDict["TraitValues"]
    traitCensCodes = dataDict["CensCodes"]
    traitCensBounds = dataDict["CensBounds"]
    patientAlleleVecs = dataDict["AlleleVecs"]
    hlaAlleles = dataDict["Alleles"]
    loci = sorted(hlaAlleles.keys())

    if hlaFileName is not None:
        hlaCountDict = fetcher.importAlleleFrequencyData(hlaFileName)
    else: ## empty dictionary
        hlaCountDict = {}

    ## choose a working folder
    work_folder = os.getcwd()
    ## get the right output folder
    if "outputFolder" in parDict.keys():
        outputFolder = os.path.join(work_folder, parDict["outputFolder"])
    else:
        outputFolder = os.path.join(work_folder, "data")
    ## TODO: make sure this folder and a figures subfolder exists

    ## get indices of forbidden alleles
    leaveOutAlleleVec = {locus : [hla in leaveOutAlleles for hla in hlaAlleles[locus]]
                         for locus in loci}

    ## find patients that should be left out
    def hasForbiddenAllele(av, locus):
        return any(aux.flatten([[a1 and a2 for a1, a2 in zip(av[i], leaveOutAlleleVec[locus])]
                                for i in range(2)]))

    leaveOutPatientIdxs = sorted(aux.flatten([[i for i, av in enumerate(patientAlleleVecs[locus])
                                              if hasForbiddenAllele(av, locus)] for locus in loci]))

    ## set the VL censoring of the left-out patients to missing
    traitCensCodes = [code if i not in leaveOutPatientIdxs else defn.missing_code
                   for i, code in enumerate(traitCensCodes)]

    ## make an ETE tree
    colorfun = lambda hlaStr: defn.locusColorDict[mhctools.MhcObject(hlaStr).locus]
    tree, nodeNames = colortrees.mkEteHlaTree(hlaNewickString, colorfun=colorfun)

    ## prepare data for Stan
    Ploidy = 2
    NumSubjects = len(traitValues)
    NumLoci = 3
    NumAlleles = [len(hlaAlleles[locus]) for locus in loci]
    AdmAlleles = mkAdmissibilityLists(patientAlleleVecs)
    NumAdmAlleles = [[[len(x) for x in xs] for xs in xss] for xss in AdmAlleles]
    NumNodes = len(nodeNames)
    treematdict = mkNodeMatrix(tree, nodeNames, hlaAlleles, represDict=represDict)
    TreeMatrix = np.concatenate([treematdict[locus] for locus in loci])
    LengthEdgeToParent = mkLengthVector(tree, nodeNames) ## specify lenfun to transform distances
    AddlAlleleData = [[hlaCountDict[hla] if hla in list(hlaCountDict.keys()) else 0
                      for hla in hlaAlleles[locus]] for locus in loci]
    ## the Stan model uses one vector for traitValues and traitCensBounds
    lrcenscodes = [defn.left_censored_code, defn.right_censored_code]
    TraitValue = [traitValues[i] if traitCensCodes[i] not in lrcenscodes else traitCensBounds[i]
                  for i in range(NumSubjects)]
    ## translate traitCensCodes to Stan censor codes
    TraitCensorType = list(map(getStanCensCode, traitCensCodes))

    ## define the multiplicity
    if multiplicity == "single":
        raise Exception("model not implemented (FIXME)") ## FIXME
    elif multiplicity == "double":
        pass ## FIXME
    else:
        raise Exception("invalid allele multiplicity")

    ## define the prior for the NodeWeights
    if prior == "norm":
        NodeWeightPrior = stantools.normal_prior_code
    elif prior == "dexp":
        NodeWeightPrior = stantools.laplace_prior_code
    else:
        raise Exception("invalid prior distribution")

    ## make data dictionary and select parameters to monitor
    data = {
        "NumLoci" : NumLoci,
        "NumAlleles" : NumAlleles,
        "AddlAlleleData" : aux.flatten(AddlAlleleData),
        "Ploidy" : Ploidy,
        "NumSubjects" : NumSubjects,
        "TraitValue" : TraitValue,
        "TraitCensorType" : TraitCensorType,
        "NumAdmAlleles" : NumAdmAlleles,
        "AdmAlleles" : aux.flatten_recursive(AdmAlleles),
        "NumNodes" : NumNodes,
        "TreeMatrix" : TreeMatrix,
        "LengthEdgeToParent" : LengthEdgeToParent,
        "WBIC" : 0,
        "NodeWeightPrior" : NodeWeightPrior
    }
    monitor = [
        "rescaledNodeWeights",
        "rescaledAlleleWeights",
        "sigmaTraitValue",
        "sigmaNodeWeight",
        "alleleFreqs",
        "admAlleleProbs",
        "traitValueLoglikes"
    ]
    ## add a specifier to the model name to indicate left out alleles
    if len(leaveOutAlleles) > 0:
        modelName += "." + ".".join(hla.Newick_str() for hla in leaveOutAlleles)
    ## choose the right jags model
    file_name = os.path.join(defn.ROOT_DIR, "stan/hla-tree-model.stan")
    ## choose a folder for the compiled stan model
    path = os.path.join(defn.ROOT_DIR, "stan/cache") ## FIXME: test if writa access (on setup: pre-compile stan models)
    ## make a Stan model object
    with open(file_name) as f:
        model_code = f.read()
    sm = stantools.CachedStanModel(model_code=model_code, path=path)
    ## run the Stan model
    if not dry_run:
        fit = sm.sampling(data=data, pars=monitor, iter=2*chain_len,
                          warmup=chain_len, thin=chain_thin, chains=num_chains)
        if verbose:
            print(fit)
        if wbic_sampling:
            data["WBIC"] = 1 ## change WBIC option
            monitor = ["sumTraitValueLoglikes"] ## only have to monitor the log likelihood
            fit_wbic = sm.sampling(data=data, pars=monitor, iter=2*chain_len,
                                   warmup=chain_len, thin=chain_thin, chains=num_chains)
            if verbose:
                print(fit_wbic)
        else:
            fit_wbic = None
    else:
        fit = None
        fit_wbic = None
    ## compute statistics, color the tree
    if fit is not None:
        chain = fit.extract()
        ## compute statistics
        WAIC = statistics.calcWAIC(aux.transpose(chain["traitValueLoglikes"]))
        if fit_wbic is not None:
            chain_wbic = fit_wbic.extract()
            WBIC = statistics.calcWBIC(chain_wbic["sumTraitValueLoglikes"])
            ## TODO: compute Rhat
        else:
            WBIC = None
        ## add the estimates to the tree
        nodeWeights = [np.array(x) for x in aux.transpose(chain["rescaledNodeWeights"])]
        colortrees.addNodeFeatures(tree, nodeNames, "beta", nodeWeights)
        colortrees.addCumulativeFeature(tree, "beta", scale_dist=True)
        expectedHlaCounts = getExpectedHlaCountsStan(chain, hlaAlleles, leaveOutPatientIdxs)
        colortrees.addLeafBars(tree, hlaAlleles, expectedHlaCounts, represDict=represDict)
    else:
        if verbose:
            print("no sample generated or no cached sample found.")
        chain = None
        WAIC = None
        WBIC = None
    ## color and plot the tree
    ## draw tree with colors indicating weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder,
      "figures/colored-mhc-tree.{}.png".format(modelName))
    colortrees.colorEteTree(tree, cfname="beta", wfname="beta", cffun=np.mean,
                            wffun=aux.getMass, cf_scale_dist=True)
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder,
      "figures/colored-mhc-utree.{}.png".format(modelName))
    if render_trees:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## draw tree with colors indicating cumulative weights
    ts = colortrees.getTreeStyle()
    figFileName = os.path.join(outputFolder,
      "figures", "colored-mhc-tree.{modelName}.cumulative.png")
    colortrees.colorEteTree(tree, cfname="beta_cumulative", wfname="beta",
                            cffun=np.mean, wffun=(lambda x: abs(np.mean(x))),
                            cf_scale_dist=False, wf_scale_dist=True)
    ## TODO legend
    ## draw a NetworkX tree
    G = colortrees.eteToNetworkX(tree)
    figFileNameNX = os.path.join(outputFolder, "figures",
        f"colored-mhc-utree.{modelName}.cumulative.png")
    if render_trees:
        if aux.isXServerAvailable():
            tree.render(figFileName, w=200, units="mm", tree_style=ts, dpi=500)
            colortrees.renderTreeNetworkX(G, filename=figFileNameNX, colorfun=colorfun)
        elif verbose:
            warnings.warn("no X-server avaliable. Tree not rendered.")
    ## return results
    rd = {
        "WAIC" : WAIC,
        "WBIC" : WBIC,
        #"tree" : tree, ## FIXME: ete3 tree can not be pickled
        "chain" : chain,
        "leaveOutPatientIdxs" : leaveOutPatientIdxs
    }
    return rd
