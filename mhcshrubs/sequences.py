from __future__ import print_function
from builtins import zip
import numpy as np
from mhcshrubs import auxiliary as aux

## data

ambiguous_aas = ['X', 'B', 'Z', 'J', '*']

## functions

def getFastaKeys(fastafilename):
    """
    just get the keys from a fasta file (TODO: obsolete?)
    """
    with open(fastafilename, 'r') as fileHandle:
        fileContent = fileHandle.read().split('>')
    ## ignore everything before the first '>'
    xss = [x.split('\n') for x in fileContent[1:]]
    keys = [xs[0].strip() for xs in xss if len(x) > 0]
    return keys[1:]


def readFastaFile(fastafilename):
    """
    import the content of a fasta file and return a dictionary (key, sequence).
    Any duplicates of a key are ignored (todo: handle duplicates in a better way)
    """
    with open(fastafilename, 'r') as fileHandle:
        fileContent = fileHandle.read().split('>')
    ## ignore everything before the first '>'
    xss = [x.split('\n') for x in fileContent[1:]]
    keys = [xs[0].strip() for xs in xss if len(xs) > 0]
    values = ["".join(xs[1:]) for xs in xss if len(xs) > 0]
    seqDict = {}
    for key, value in zip(keys, values):
        if key not in list(seqDict.keys()):
            seqDict[key] = value
    return seqDict


def mkPseudoSequence(fullSeq, subSeqPos, start=0):
    """
    make a 'psuedo' sequence (i.e. a non-contiguous sub-sequence of a full sequence.
    input:
    - fullSeq -- the full source sequence
    - subSeqPos -- the positions in the full sequence belonging to the sub-sequence
    - start -- subSeqPos starts counting at the start position of fullSeq
    output: the pseudo-sequence
    """
    return "".join([fullSeq[i+start] for i in subSeqPos if i < len(fullSeq)])


def calcSequenceDist(s1, s2, aaDistDict=None):
    """
    compute the distance between two aa sequences, using a vectorspace representation
    of the amino-acids.
    input:
    - s1 -- the first sequence
    - s2 -- the second sequence
    - aaDistDict -- a dictionary string -> (string -> float)
      If aaDistDist==None, the Hamming distance is computed (default: None)
    output:
    - dist -- the calculated distance
    TODO:
    - handle ambiguous aa's
    - pass a variable for different norms (Euclidean, ell_p etc.)
    """
    assert len(s1) == len(s2), "Sequences must have the same length"
    if aaDistDict is None:
        ## count the number of non-matching aa's
        return aux.countMatching(lambda a_b: a_b[0]!=a_b[1], list(zip(s1, s2)))
    else:
        ## use the aaDistDict
        return np.sqrt(sum(aaDistDict[a1][a2]**2 for a1, a2 in zip(s1, s2)))


def mkAaDistDict(aaCovFileName, comment_char='#', use_aas=None, ignore_aas=None,
                 center_covmat=False, rescale=False):
    """
    make a amino-acid distance matrix (represented by a dictionary)
    using the contents of the given file.
    The distance matrix (dictionary) can be used in e.g. calcSequenceDist.
    The file must contain a table with the following format
      A C D E ...
    A 1 2 4 2 ...
    C 2 1 2 3 ...
    : : : : :
    with any whitespace as separator.
    input:
    - aaCovFileName -- the name of the file containing the matrix
    - comment_char -- lines starting with this character are ignored (default: #)
    - use_aas -- use these 1-character amino-acid codes.
      if use_aas=None, the characters in the file are considered (default: None)
    - ignore_aas -- ignore these 1-character amino-acid codes.
      if ignore_aas=None, all characters in the file are used.
      ignore_aas overrules use_aas (default: None).
    - center_covmat -- apply a transformation that 'centers' the covariance matrix.
      C -> H^t C H with H = I - J/n, I is the identify, J is a block of ones,
      and n is the dimension.
      TODO: is this trivial? computing distance appears to undo centering
    - rescale -- rescale the distances such that the average distance is 1 (default: False)
    output:
    - distDict -- a dictionary of dictionaries with distances
    """
    with open(aaCovFileName, 'r') as f:
        table = [row.split() for row in f.read().split('\n') if row != '' and row[0] != comment_char]
    aas = table[0]
    if use_aas is not None:
        if any(a not in aas for a in use_aas):
            raise Exception("required amino acid is not listed in the covariance matrix file")
        aas = use_aas
    if ignore_aas is not None:
        aas = [aa for aa in aas if aa not in ignore_aas]
    covDict = dict((row[0], dict((aa, float(x)) for aa, x in zip(aas, row[1:]))) for row in table[1:])
    ## convert covariances to distances
    covMat = np.array([[covDict[a1][a2] for a1 in aas] for a2 in aas])
    if center_covmat:
        n = len(aas)
        H = np.eye(n) - np.ones((n,n))/n
        covMat = np.dot(H.T, np.dot(covMat, H))
    distMat = aux.covToDistMatrix(covMat)
    if rescale:
        distMat /= np.mean(distMat)
    distDict = dict((a1, dict((a2, distMat[i1,i2]) for i2, a2 in enumerate(aas))) for i1, a1 in enumerate(aas))
    return distDict


def mkEquivalenceClasses(sequenceDict):
    """
    based on sequence identity, make equivalence classes of keys.
    input:
    - sequenceDict -- dictionary key -> sequence
    output:
    - representatives -- a system of representatives
    - represDict -- dictionary key -> representative
    TODO: pass rules for choosing a representative?
    """
    ## find unique sequences
    useqs = aux.unique(list(sequenceDict.values()))
    ## make equivalence classes
    equivclasses = [[key for key, seq in sequenceDict.items() if seq == useq] for useq in useqs]
    ## pick the "smallest" key as a representative
    representatives = [min(equivclass) for equivclass in equivclasses]
    ## return a dict key -> representative
    represDict = dict(aux.flatten([[(key, repres) for key in equivclass]
                                   for repres, equivclass in zip(representatives, equivclasses)]))
    return (representatives, represDict)
