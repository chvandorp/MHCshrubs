from __future__ import print_function
from builtins import (next, map, filter)
import scipy.cluster.hierarchy as sclush
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import multiprocessing
import ete3
import os
import itertools
from mhcshrubs import auxiliary as aux
from mhcshrubs import definitions as defn
from mhcshrubs import (mhctools, mhcbind, sequences)


def mkPseqDict(pSeqFileName): ## MHC_pseudo.dat
    """
    input:
    - pSeqFileName -- name of file with pseudo sequences
    output:
    - pSeqDict -- dictionary with HLA alleles -> pseudo sequences
    """
    classic_human_loci = ["HLA-A", 'HLA-B', "HLA-C"]
    with open(pSeqFileName, 'r') as pSeqFileHandle:
        pSeqDict = pSeqFileHandle.read().split('\n')
    pSeqDict = [row.split(' ') for row in pSeqDict]
    pSeqDict = dict((mhctools.MhcObject(mhctools.fromNetMhcFormat(row[0])), row[1])
                    for row in pSeqDict if any([x in row[0] for x in classic_human_loci]))
    return pSeqDict


def mkPseqFastaFile(pSeqDict, hlaAlleles, outFileName):
    with open(outFileName, 'w') as fileHandle:
        for X in sorted(hlaAlleles.keys()):
            for hla in hlaAlleles[X]:
                th = hla.Newick_str()
                sh = pSeqDict[hla]
                fileHandle.write(">{0}\n{1}\n".format(th, sh))
    ## end of 'with' block closes file

def getNewick(node, newick, parentdist, leaf_names):
    """
    Transform scipy.cluster.hierarchy trees into newick format.
    Function builds the newick string recursively.
    Thanks to "jfn" (http://stackoverflow.com/users/3922554/jfn).
    input:
    - a node in the hierarchy tree
    - the newick string so-far
    - the distance to the parent
    - the names of the leafs
    output: a newick string
    """
    if node.is_leaf():
        return "{0:s}:{1:f}{2:s}".format(leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):{0:f}{1:s}".format(parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",{:s}".format(newick), node.dist, leaf_names)
        newick = "({:s}".format(newick)
        return newick

def mkDistanceMatrix(sequenceDict, aaDistDict=None):
    """
    Make a distance matrix for the sequences in the dictionary.
    input:
    - sequenceDict -- a dictionary of sequences
    - aaDistDict -- a dictionary for distances between amino-acids
      if None is passed, then the Hamming distance is used (default: None)
    output:
    - distmat -- a distance matrix
    - distmat_header -- the correct order of the sequences
      (using the keys of sequenceDict)
    """
    distmat_header = sorted(sequenceDict.keys())
    n = len(distmat_header)
    distmat = np.zeros((n,n))
    for i1, s1 in enumerate(distmat_header):
        for i2, s2 in enumerate(distmat_header):
            if i1 > i2:
                d = sequences.calcSequenceDist(sequenceDict[s1], sequenceDict[s2], aaDistDict=aaDistDict)
                distmat[i1, i2] = distmat[i2, i1] = d
    return (distmat, distmat_header)

def collapseDistanceMatrix(distmat, header, epsilon=1e-5):
    """
    remove duplicate rows and columns, make a dictionary of representatives
    """
    n, m = distmat.shape
    if n != m:
        raise Exception("can not collapse non-square matrix")
    represIdxDict = {}
    for i in range(n):
        for j in range(i,n):
            if distmat[i,j] <= epsilon:
                if j not in represIdxDict.keys():
                    represIdxDict[j] = i
    ## reduce the matrix
    represIdxs = aux.unique(represIdxDict.values())
    cdistmat = distmat[np.ix_(represIdxs, represIdxs)]
    cheader = [header[i] for i in represIdxs]
    represDict = {header[i] : header[j] for i, j in represIdxDict.items()}
    return (cdistmat, cheader, represDict)


def scaleHClus(Z, s=1.0):
    """
    Scale the depth of the hierarchical clustering Z to s

    Args:
        Z -- clustering returned by scipy.cluster.hierarchy.linkage

    Kwargs:
        s -- target depth (default = 1)

    Returns:
        scaled copy of Z
    """
    NZ = Z
    NZ[:,2] = s * Z[:,2] / np.max(Z[:,2])
    return NZ


def mkNewickTree(sequenceDict, aaDistDict=None, collapse=False):
    """
    Turn a list of sequences (dictionary key -> seq) into a newick format tree.
    The aaDistDict is used for distances between amino-acids.
    The keys in the sequence dictionary are used as leaf names
    input:
    - sequenceDict -- a dictionary of sequences
    - aaDistDict -- a dictionary for distances between amino-acids
    - collapse -- collapse equivalence classes to one representative (default: False)
    output:
    - newick_str -- a tree in newick format
    - represDict -- dictionary mapping alleles to their representatives
      that are used as leafs in the tree. trivial if collapse is False
    """
    if collapse:
        ## find equivalence classes of sequences
        representatives, represDict = sequences.mkEquivalenceClasses(sequenceDict)
        ## restrict the sequence dict to just the representatives
        sequenceDict = {key : sequenceDict[key] for key in representatives}
    else:
        represDict = {k : k for k in sequenceDict.keys()} ## trivial
    ## make a newick tree using the (representatives of the) sub-sequences and the aaDistDict
    distmat, distmat_header = mkDistanceMatrix(sequenceDict, aaDistDict=aaDistDict)
    ## flatten the distance matrix (TODO: use scipy function)
    flat_distmat = aux.mkFlatDistMat(distmat)
    ## use scipy.cluster.hierarchy to make a dendogram
    Z = sclush.linkage(flat_distmat, method='average') ## UPGMA, cf. MhcCluster
    ## rescale Z
    Z = scaleHClus(Z)
    ## TESTING: plot the dendogram
    #fig, (ax1) = plt.subplots(1, 1, figsize=(20,5))
    #sclush.dendrogram(Z, orientation='top', labels=distmat_header, ax=ax1)
    #plt.setp(ax1.get_xticklabels(), rotation=90)
    #fig.savefig("test_dendrogram.png")
    #plt.close(fig)
    ## END TESTING
    tree = sclush.to_tree(Z, False)
    ## convert the distmat_header into values that are allowed in a newick tree
    nw_distmat_header = [x.Newick_str() for x in distmat_header]
    newick_str = getNewick(tree, "", tree.dist, nw_distmat_header)
    return (newick_str, represDict)


def mkClustaloDistanceMatrix(fastaFileName, verbose=False):
    """
    sequences in a fasta file are aligned with clustalo, which also
    computes a distance matrix. This matrix is loaded from a file
    and returned, together with a header indicating the order of the rows
    and columns.
    input:
    - fastaFileName -- name of a fasta file
    output:
    - distmat -- the matrix (numpy array)
    - distmat_header -- the sequence IDs from the fasta keys (> ID)
    """
    logFileName = aux.mkUniqueFile("/tmp/clustalo", "log")
    logFileHandle = open(logFileName, 'w')
    errFileName = aux.mkUniqueFile("/tmp/clustalo", "err")
    errFileHandle = open(errFileName, 'w')
    distmatFileName = aux.mkUniqueFile("/tmp/distance_matrix", "txt")

    retcode = subprocess.call(["clustalo",
                               "--full", ## make sure distance matrix is computed
                               "--force", ## overwrite existing files
                               "--infile={}".format(fastaFileName),
                               #"--profile1={}".format(fastaFileName),
                               "--is-profile",
                               "--distmat-out={}".format(distmatFileName)],
                               stdout=logFileHandle, stderr=errFileHandle)
    logFileHandle.close()
    errFileHandle.close()

    if verbose:
        logFileHandle = open(logFileName, 'r')
        print(logFileHandle.read())
        logFileHandle.close()
        errFileHandle = open(errFileName, 'r')
        print(errFileHandle.read())
        errFileHandle.close()

    distmatFileHandle = open(distmatFileName, 'r')
    distmat = distmatFileHandle.read().split('\n')
    distmatFileHandle.close()
    distmat = [row.split() for row in distmat[1:] if row != '']
    distmat_header = [row[0] for row in distmat]
    distmat = np.array([list(map(float, row[1:])) for row in distmat])
    return (distmat, distmat_header)


def mkClustaloNewickTree(fastaFileName, verbose=False):
    """
    make a newick tree using sequences in a fasta file.
    first, use mkDistanceMatrix to compute a distance matrix between the
    sequences, then use the scipy.cluster module to compute a tree.
    This tree is then converted into a newick string with the function
    getNewick
    input:
    - fastaFileName -- file with sequences
    - verbose -- print messages
    output:
    - newick_str -- the newick string representing the tree
    """
    ## get the distance matrix
    distmat, distmat_header = mkClustaloDistanceMatrix(fastaFileName, verbose=verbose)
    flat_distmat = aux.mkFlatDistMat(distmat)
    ## use scipy.cluster.hierarchy to make a dendogram
    Z = sclush.linkage(flat_distmat, method='average') ## UPGMA, cf. MhcCluster
    ## scale Z such that depth is 1
    Z = scaleHClus(Z)
    ## optionally plot the dendogram
    #fig, (ax1) = plt.subplots(1, 1, figsize=(20,5))
    #sclush.dendrogram(Z, orientation='top', labels=distmat_header, ax=ax1)
    #plt.setp(ax1.get_xticklabels(), rotation=90)
    tree = sclush.to_tree(Z, False)
    newick_str = getNewick(tree, "", tree.dist, distmat_header)
    return newick_str


def mkTrivialNewickTree(hlaAlleles, incl_null_alleles=False):
    """
    Make a completely star-shaped tree in newick format

    Args:
        hlaAlleles -- a dict of locus -> list of MhcObject-s

    Kwargs:
        incl_null_alleles -- if False (default), ignore Null alleles

    Returns:
        newick_str -- a string representing a tree in the newick format,
            with each allele connected to the root.

    Todo:
        allow for equivalence classes
    """
    hlastrs = [hla.Newick_str() for X in sorted(hlaAlleles.keys())
               for hla in hlaAlleles[X] if hla.suffix != 'N' or incl_null_alleles]
    newick_str = "(" + ",".join("{}:1.0".format(h) for h in hlastrs) + ");"
    return newick_str


def mkGroupNewickTree(hlaAlleles, incl_null_alleles=False, starlikeness=0.5):
    """
    Make a newick tree representing one-field level typing

    Args:
        hlaAlleles -- a dict of locus -> list of MhcObject-s

    Kwargs:
        incl_null_alleles -- if False (default), ignore Null alleles
        starlikeness -- if 1.0, the group (type) branches have length 0.0,
            and hence the tree is star-shaped (of depth 1.0).
            If 0.0, all subtype branches have length 0.0.
            The tree is again star-shaped, but the leafs represent
            clusters of subtypes (i.e. the types). Default: 0.5

    Returns:
        newick_str -- a newick tree with intermediate nodes at the one-field level
    """
    if starlikeness > 1.0 or starlikeness < 0.0:
        raise ValueError("starlikeness should be between 0 and 1")
    alleles = [hla for X in sorted(hlaAlleles.keys())
               for hla in hlaAlleles[X] if hla.suffix != 'N' or incl_null_alleles]
    groups = aux.unique(hla.type for hla in alleles)
    clusters = [[hla for hla in alleles if hla in grp] for grp in groups]
    ## sub-trees for each group
    cluster_strs = ["(" + ",".join("{0:s}:{1:f}".format(h.Newick_str(), starlikeness)
                    for h in clus) + ")"
                    for clus in clusters]
    ## combine sub-trees in a group-based tree
    newick_str = "(" + ",".join("{0:s}:{1:f}".format(n, 1.0-starlikeness)
                                for n in cluster_strs) + ");"
    return newick_str


def treeToCovmat(newickTree):
    """
    Use a newick tree to produce a covariance matrix.

    The covariance between alleles is defined as the distance from root
    to the MRCA of the two allele leafs.

    Args:
        newickTree -- the newick tree with HLA names at the leafs (tree format)

    Returns:
        a tuple consisting of:
            - covmat -- the resulting covariance matrix
            - names -- a list of leaf names in the used order
    """
    tree = ete3.Tree(newickTree)
    root = tree.get_tree_root()
    ## get all alleles in the tree
    names = sorted([leaf.name for leaf in tree])
    H = len(names)
    covmat = np.zeros((H, H))
    for leaf1 in tree:
        idx1 = names.index(leaf1.name)
        for leaf2 in tree:
            idx2 = names.index(leaf2.name)
            ca = tree.get_common_ancestor(leaf1, leaf2)
            covmat[idx1, idx2] = ca.get_distance(root) ## TODO: should this be the square?
    return (covmat, names)



## methods for peptide-binding based clustering

def mkPeptideBindingNewickTree(fastaFileName, hlaAlleles, collapse=False,
                               verbose=False, dry_run=False, use_cache=False,
                               cache_dir="/tmp/", parallel=True):
    """
    Make a newick tree encoding HLA-similarities based on peptide binding.

    Args:
        fastaFileName -- proteines of the pathogen of interest
        hlaAlleles -- alleles of interest (sorted by locus)

    Kwargs:
        collapse -- if True, choose a representative for identical alleles
        verbose -- print messages? (default: False)
        dry_run -- don't run NetMHCpan (default: False)
        use_cache -- don't run NetMHCpan, but use cached data
        cache_dir -- directory for caching NetMHCpan output
        parallel -- if True, use multiple cpu cores

    Returns:
        a tuple containing:
            - the tree in newick format
            - a dictinary for allele representatives

    """
    mhcList = [hla for X in sorted(hlaAlleles.keys()) for hla in hlaAlleles[X]]
    netMHCpanPath = defn.NETMHCPAN_PATH

    ## determine the number of cpu cores to use
    if parallel:
        max_workers = max(1, multiprocessing.cpu_count()-1)
    else:
        max_workers = 1
    ## call NetMHCpan
    dry_run_netmhc = dry_run or use_cache
    results = mhcbind.callNetMHCpanParallel(fastaFileName, mhcList, netMHCpanPath,
                                            cache_dir=cache_dir, threads=max_workers,
                                            verbose=verbose, dry_run=dry_run_netmhc)
    ## collect filenames
    predFileNames = [x[1] for x in results]
    ## collect the output of NetMHCpan
    predDict = {}
    for predFileName in predFileNames:
        predDict.update(mhcbind.importNetMhcOutput(predFileName))
    ## make a correlation matrix
    corrMat, mhcList = mhcbind.mkBindingCorrelationMatrix(predDict, var='rank', fun=np.log)
    ## convert the correlation matrix into a distance matrix
    distMat = aux.covToDistMatrix(corrMat)
    ## collapse the distance matrix and choose representatives
    if collapse:
        distMat, mhcList, represDict = collapseDistanceMatrix(distMat, mhcList)
    else:
        represDict = {mhc : mhc for mhc in mhcList} ## trivial
    ## and make a flat version for scipy.cluster
    flatDistMat = aux.mkFlatDistMat(distMat)
    ## cluster the alleles and transform to newick tree format
    mhc_names = [mhc.Newick_str() for mhc in mhcList]
    Z = sclush.linkage(flatDistMat, method="average") ## UPGMA, cf. MhcCluster
    ## normalize Z to depth 1
    Z = scaleHClus(Z)
    tree = sclush.to_tree(Z, False)
    newick_str = getNewick(tree, "", tree.dist, mhc_names)
    return (newick_str, represDict)


## make a tree based on specific parts of the MHC sequence (e.g. the KIR epitope)

def mkSubSeqNewickTree(hlaAlleles, positions, start=0, aaDistDict=None, verbose=False):
    """
    make a tree using particular positions in the MHC sequence.
    example: re-create netMHCpan pseudo sequences.
    input:
    - hlaAlleles -- the alleles at the leafs of the tree
    - fastaFileName -- aligned MHC protein sequences
    - positions -- the positions of interest
    - start -- the alignment could start at a different position
      (see the function sequences.mkPseudoSequence)
    - verbose -- print comments (default: False) TODO: unused
    output:
    - newick_str -- the resulting tree in Newick format
    - represDict -- a dictionary mapping the alleles to the chosen representative
    """
    ## import IMGT sequences
    imgtSeqDict = {}
    work_folder = os.getcwd() ## FIXME pass path and name of fasta files?
    loci = sorted(hlaAlleles.keys())
    for locus in loci:
        imgtFileName = os.path.join(work_folder, "data/hla/HLA-{}-IMGT-protein-alignment.fasta".format(locus))
        imgtSeqDict.update(sequences.readFastaFile(imgtFileName))
    ## get sub-sequences
    subSeqDict = {}
    imgtKeys = sorted(imgtSeqDict.keys())
    for locus in loci:
        for hla in hlaAlleles[locus]:
            ## find a key in imgtSeqDict containing the HLA allele
            key = next(filter(lambda k: mhctools.MhcObject(k) in hla, imgtKeys), None)
            ## NB: changed behavior: previously hla.subtype_str() in key
            if key is None:
                raise Exception("HLA allele {} not present in IMGT database".format(hla))
            seq = imgtSeqDict[key]
            subSeq = sequences.mkPseudoSequence(seq, positions, start)
            subSeqDict[hla] = subSeq
    ## make a newick tree, collapsing equivalence classes into one representative in the tree
    newick_str, represDict = mkNewickTree(subSeqDict, aaDistDict=aaDistDict, collapse=True)
    if verbose:
        print("representative alleles:")
        for repres in aux.unique(list(represDict.values())):
            print(repres, subSeqDict[repres])
    return (newick_str, represDict)
