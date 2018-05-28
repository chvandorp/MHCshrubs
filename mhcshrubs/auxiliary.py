from __future__ import (print_function, division)
from builtins import (map, range)
import numpy as np
import scipy.stats as sts
import string
import os
import sys
from subprocess import Popen, PIPE ## for X11 check

missing_code = "missing" ## TODO: parameters
left_censored_code = "left_censored" ## the given value is an upper bound
right_censored_code = "right_censored" ## the given value is a lower bound
uncensored_code = "uncensored"

auxiliaryLowerCensBound = 0.0 ## TODO: used in JAGS. is this nescesary?

def fst(pair): (a, _) = pair; return a

def snd(pair): (_, b) = pair; return b

def transpose(xss):
    n = len(xss[0])
    return [[xs[i] for xs in xss] for i in range(n)]

def invertRepresDict(represDict):
    """
    Turn a dictionary of representatives into equivalence classes,
    indexed by the representatives.
    input:
    - represDict -- a dictionary with representatives
    output:
    - equivClassDict -- a dictionary with equivalence classes
    """
    representatives = set(represDict.values())
    equivClassDict = dict((rep, [key for key, val in list(represDict.items()) if val == rep])
                          for rep in representatives)
    return equivClassDict



def flatten(xss): return [x for xs in xss for x in xs]

def flatten_recursive(xs):
    if type(xs) != list:
        return [xs]
    ## else: we have a list...
    fxs = []
    for x in xs:
        fxs += flatten_recursive(x)
    return fxs

def un_flatten(Nums, FlatData):
    """
    convert a flat data structure (a list) back to a list-of-lists
    TODO: allow for a flexible depth (now always 3)
    """
    SortedData = []
    i = 0
    for s, nss in enumerate(Nums):
        xss = []
        for p, ns in enumerate(nss):
            xs = []
            for ell, n in enumerate(ns):
                x = FlatData[i:i+n]
                i += n
                xs.append(x)
            xss.append(xs)
        SortedData.append(xss)
    return SortedData


def unique(xs):
    """Returns sorted, unique elements of xs"""
    return sorted(list(set(xs)))

def split_list(xs, n):
    """
    Split the list xs into n parts of (approx) equal length.
    The number n is a maximum: when len(xs) < n, split_list returns
    [[x] for x in xs] instead.
    when len(xs) == n*q + r with r >= 0 and q > 0 maximal,
    the first r chunks are of length q+1 and the remaining chunks of length q
    """
    r = len(xs) % n
    q = len(xs) // n ## integer division
    if r == 0:
        return [xs[q*i:q*(i+1)] for i in range(n)]
    ## get r chunks of size q+1 and n-r chunks of size q
    return split_list(xs[:r*(q+1)], r) + split_list(xs[r*(q+1):], n-r)


def countMatching(condition, seq):
    """Returns the amount of items in seq that return true from condition"""
    return sum(1 for item in seq if condition(item))


def basicTsvImport(fileName, sep='\t', unquote=False, has_header=False):
    fileHandle = open(fileName, 'r')
    table = fileHandle.read().split('\n')
    fileHandle.close()
    table = [row.split(sep) for row in table if row != '']
    if unquote:
        def unquoter(string): return [ch for ch in string if ch != '"']
        table = [list(map(unquoter, row)) for row in table]
    if has_header:
        header = table[0]
        table = table[1:]
        return (table, header)
    else:
        return table

def getMass(xs): ## mass above (or below) zero
    massAboveZero = len([x for x in xs if x > 0]) / len(xs) ## division results in float
    if np.mean(xs) > 0: return massAboveZero
    else: return 1-massAboveZero

def log_sum_exp(xs):
    if len(xs) == 0:
        return np.nan
    ## else, we can take the max...
    m = np.max(xs)
    return m + np.log(np.sum([np.exp(x-m) for x in xs]))

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def is_semi_pos_def(x):
    return np.all(np.linalg.eigvals(x) >= 0)


def is_distance_matrix(D, verbose=False):
    ## check triangle inequality
    ## TODO: see also scipy functions as scipy.spatial.distance.is_valid_dm
    n, m = D.shape
    ok = True
    if n != m: return False
    for i in range(n):
        ## distance to self is 0?
        DistSelfOk = (D[i,i] == 0)
        if not DistSelfOk:
            if verbose: print("non zero diagonal", i, ":", D[i,i])
            ok = False
        for j in range(n):
            ## positive?
            PosOk = (D[i,j] >= 0)
            if not PosOk:
                if verbose: print("not positive", i, j, ":", D[i, j])
                ok = False
            ## symmetric?
            SymOk = (D[i,j] == D[j,i])
            if not SymOk:
                if verbose: print("not symmetric", i, j, ":", D[i, j], D[j, i])
                ok = False
            for k in range(n):
                ## triangle inequality
                TrEqOk = ((D[i,j] + D[j,k]) >= D[i,k])
                if not TrEqOk:
                    if verbose: print("invalid triangle inequality", i, j, k, ":", D[i,j], D[j, k], D[i, k])
                    ok = False
    return ok

def mkFlatDistMat(distmat):
    n = len(distmat)
    M = np.array(distmat)
    return M[np.triu_indices(n, k=1)]


def mkRDmatrix(D):
    """
    Auxiliary function for converting distance matrices into points in Euclidean space.

    This can be used to compute correlation matrices.
    Suppose that we have points x_i in R^n, then the distance matrix D
    satisfies D_ij = \|x_i - x_j\|. We now define M such that
    M_ij = <x_i - x_1 | x_j - x_1> = <y_i | y_j> with y_i = x_i - x_1.
    We have M = Y Y^t with Y = (...,y_i,...).
    The matrix M has the form M_ij = 1/2*(D_1j^2 + D_i1^2 - D_ij^2),
    because the latter equals:
    1/2*(\|x_1 - x_j\|^2 + \|x_i - x_1\|^2 - \|x_i - x_j\|^2) =
    1/2*(\|x_1\|^2 + \|x_j\|^2 - 2<x_1|x_j>
       + \|x_i\|^2 + \|x_1\|^2 - 2<x_i|x_1>
       - \|x_i\|^2 - \|x_j\|^2 + 2<x_i|x_j>) =
    <x_i|x_j> - <x_i|x_1> - <x_1|x_j> + <x_1|x_1> =
    <x_i - x_1 | x_j - x_1>    QED.

    Args:
        D -- distance matrix

    Returns:
        given a distance matrix D = D_ij,
        return M_ij = 1/2*(D_1j^2 + D_i1^2 - D_ij^2)
    """
    n, m = D.shape
    if n != m:
        raise Exception("Matrix is not square")
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            M[i, j] = 0.5*(D[i,0]**2 + D[0,j]**2 - D[i,j]**2)
    return M

def mkGramianMatrix(M):
    """
    find a matrix X such that M = X^tX. M must be semi positive definite
    in which case we can use the eigen decomposition M = U S U^t.
    Then we get X = sqrt(S)U^t
    """
    if not is_semi_pos_def(M):
        raise Exception("Matrix is not semi positive-definite")
    evals, U = np.linalg.eigh(M)
    S = np.diag(np.sqrt(evals))
    X = np.dot(S, U.transpose())
    return X

def mkCorrMatrix(hlaDistMat, hlaDistMatHeader, hlaAlleles):
    """
    Make a correlation matrix from a distance matrix.

    This is a bit problematic, because there is no obvious way to do this
    However, assuming that data X is normalized and centralized, we have
    Cor(X_i,X_j) = Cov(X_i, X_j) = E(X_iX_j) = 1/n <X_j|X_j>, and for the
    Euclidean distance, we get
    \|X_i - X_j\|^2 = \|X_i\|^2 + \|X_j\|^2 - 2 <X_i|X_j> = 2(1 - <X_i|X_j>)
    Therefore, we have Cor(X_i, X_j) = 1 - 1/2 * Dist(X_i, X_j)^2.
    Whenever the original data was NOT normalized, we can take the following
    approach: we reduce the n x n distance matrix D to RD using mkRDmatrix,
    then, the RD matrix is decomposed into X with mkGramianMatrix.
    The correlations are computed from with X

    Args:
        distmat -- a given distance matrix

    Returns:
        corrmat -- correlation matrix
    """
    ## transform distance matrix to covariance matrix
    if not aux.is_distance_matrix(hlaDistMat, verbose=False):
        raise Exception("Ill-defined distance matrix")
    M = aux.mkRDmatrix(hlaDistMat)
    X = aux.mkGramianMatrix(M)
    corrmat = np.corrcoef(X.transpose())
    return corrmat

def covToDistMatrix(covMat):
    """
    \|x-y\|^2 = (n-1) * (<x|x> + <y|y> -2<x|y>)
    notice that Cov(x, y) = 1/(n-1) * <x - <x>|y - <y>>
    """
    n, m = np.shape(covMat)
    if n != m:
        raise Exception("covariance matrix is not square")
    distMat = [[np.sqrt((n-1)*(covMat[i,i] + covMat[j,j] - 2*covMat[i,j]))
                for i in range(n)] for j in range(m)]
    return np.array(distMat)

def getRanStr(n):
    """
    returns a random string of length n.
    can be used for random file name formation
    """
    characters = string.ascii_uppercase + string.digits + string.ascii_lowercase
    return "".join(np.random.choice(list(characters), size=n))

def mkUniqueFile(fileBase, fileExt, n=10):
    """
    given a file-base and an extension, create a unique file name
    by adding some (n) random characters.
    If by a miracle, the file already exists, a new random string is tried.
    The file name is returned
    """
    createdNewFile = False
    while not createdNewFile:
        fileName = "{0}-{1}.{2}".format(fileBase, getRanStr(n), fileExt)
        try:
            fileHandle = open(fileName, 'x')
            fileHandle.close() ## we only need to "touch" the file
            createdNewFile = True
        except IOError:
            ## we expect an error if the file exists.
            pass
    return fileName

def isXServerAvailable():
    """
    Test if X11 is available.

    ETE needs an X server for rendering trees. When running the scripts
    remotely, an X server might not be available, and in this case ETE
    will crash Python. Use this function to test if an X server is available
    before rendering trees etc. On servers: run the scripts with xvfb-run

    Returns:
        boolean
    """
    dis = os.environ.get("DISPLAY", None)
    if dis is not None:
        p = Popen(["xset", "-q"], stdout=PIPE, stderr=PIPE)
        p.communicate()
        return p.returncode == 0
    else:
        return False

def posPoisson(mu):
    """
    Return a Poisson-distributed random variable, conditioned on being > 0.

    Args:
        mu -- the expectation

    Returns:
        a positive integer
    """
    x = 0
    while x == 0:
        x = sts.poisson.rvs(mu)
    return x

def truncNorm(mu, sigma, left=None, right=None):
    """
    A truncated normal distribution.

    Similar to scipy.stats version,
    but here the boundaries do not require transformation,
    and the user can choose which side(s) is (are) truncated.

    Args:
        mu -- mean of the non-truncated distribution
        sigma -- standard deviation of the non-truncated distribution

    Kwargs:
        left -- (optional) left boundary. If None, left = -infinity (default: None)
        right -- (optional) right boundary. If None, right = infinity (default: None)

    Returns:
        random deviate
    """
    if left is not None and right is not None:
        if left == right:
            return right
        elif left > right:
            raise Exception("left boundary is larger than right boundary")
    while True:
        x = sts.norm.rvs(loc=mu, scale=sigma)
        if (left is None or x >= left) and (right is None or x <= right):
            return x


def randomGroupMasking(n, p=0.5):
    """
    Make a partition of (1, 2, ..., n).

    Groups are labeled 1,2,..,g (TODO: return group label)
    the parameter p determines the bias towards
    many small groups (p = 1) or few large groups (p = 0)

    Args:
        n -- number of items

    Kwargs:
        p -- small group bias

    Returns:
        dict of item to group
    """
    xs = np.arange(1,n+1)
    np.random.shuffle(xs) ## in place shuffle
    grps = []
    grp = []
    for x in xs:
        grp.append(x)
        if sts.bernoulli.rvs(p) == 1:
            grps.append(grp)
            grp = []
    ## after the for loop, there may still be a new group to add...
    if len(grp) > 0:
        grps.append(grp)
    ## map integers to groups
    itog = {}
    for g, grp in enumerate(grps):
        for i in grp:
            itog[i] = grp
    return itog

def prob_any(ps):
    return 1.0 - np.prod(1.0 - np.array(ps))
