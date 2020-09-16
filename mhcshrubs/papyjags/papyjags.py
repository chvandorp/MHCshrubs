"""@package papyjags
This is Possibly Another PYthon interface for JAGS.

The object JagsModel was designed with pystan.StanModel in mind.

Examples:
    >>> jm = JagsModel(model_code="model {X ~ dnorm(0,1)}")
    defines a Standard Normal random variable X. The statement
    >>> chain, = jm.sampling(iter=100, pars=["X"])
    produces a 100 samples of X. They can be retrieved by
    >>> xs = chain["X"]
"""
from __future__ import (print_function, division)
from builtins import (zip, str, map, range, object)
import subprocess ## for calling jags
import time ## for sleep
import numpy as np
import numpy.random as rand
import scipy.stats as sts
import matplotlib.pyplot as plt
import string ## for random identifiers
import os ## for creating directories
import os.path ## for joining filenames
import glob ## for checking existing files
import sys
import tqdm ## for printing progressbars
import re ## for parsing JAGS output
import warnings
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

## aux functions
def flatten(xss): return [x for xs in xss for x in xs]
def fst(pair): (a, _) = pair; return a
def snd(pair): (_, b) = pair; return b
def transpose(xss): return [[xs[i] for xs in xss] for i in range(len(xss[0]))]

def isJagsInstalled():
    """Check whether `jags` is in the PATH and marked as executable"""
    from shutil import which
    return which("jags") is not None

## functions to simulate the R-dump-file
def RdumpScalar(name, x):
    x_str = str(x) if not np.isnan(x) else 'NA'
    string = "`{0}` <- {1}\n".format(name, x_str)
    return string

def RdumpVector(name, xs):
    xs_str = [str(x) if not np.isnan(x) else "NA" for x in xs]
    string = "`{0}` <- c({1})\n".format(name, ", ".join(xs_str))
    return string

def RdumpMatrix(name, xss):
    n = len(xss)
    m = len(xss[0])
    xs = flatten(xss)
    xs_str = [str(x) if not np.isnan(x) else "NA" for x in xs]
    template = "`{0}` <- structure(c({1}), .Dim=c({2}, {3}))\n"
    string = template.format(name, ", ".join(xs_str), m, n)
    return string

## functions to create the necessary files for JAGS

def mkJagsDataFile(path, fileBaseName, pd, override):
    datFileName = os.path.join(path, "{0}.dat".format(fileBaseName))
    if not override and os.path.isfile(datFileName):
        return datFileName
    ## else, create the file
    with open(datFileName, 'w') as datFileHandle:
        for p, d in pd.items():
            data = np.array(d)
            if len(data.shape) == 0: ## data is a scalar
                datFileHandle.write(RdumpScalar(p, data))
            elif len(data.shape) == 1: ## data is a vector
                datFileHandle.write(RdumpVector(p, data))
            elif len(data.shape) == 2: ## data is a matrix
                datFileHandle.write(RdumpMatrix(p, data))
            else:
                raise Exception("n-dim array with n>2 not implemented")
    return datFileName

def mkJagsParameterFile(path, fileBaseName, pd, chain, run, override):
    parFileName = os.path.join(path, "{0}.run{1}.chain{2}.par".format(fileBaseName, run, chain))
    if not override and os.path.isfile(parFileName):
        return parFileName
    ## else, create the parameter file
    with open(parFileName, 'w') as parFileHandle:
        ## set a seed for the random number generator
        seed = rand.randint(1,1000000)
        parFileHandle.write("`.RNG.name` <- \"base::Mersenne-Twister\"\n")
        parFileHandle.write("`.RNG.seed` <- {:d}\n".format(seed))
    return parFileName


def mkJagsScriptFile(path, fileBaseName, pars, iter, warmup, thin, chains, run, override):
    jgsFileName = os.path.join(path, "{0}.run{1}.jgs".format(fileBaseName, run))
    if not override and os.path.isfile(jgsFileName):
        return jgsFileName
    ## else, create the file
    with open(jgsFileName, 'w') as jgsFileHandle:
        jgsFileHandle.write("cd {}\n".format(path))
        jgsFileHandle.write("model in {}.bug\n".format(fileBaseName))
        jgsFileHandle.write("data in {}.dat\n".format(fileBaseName))
        jgsFileHandle.write("compile, nchains({})\n".format(chains))
        for chain in range(1, chains+1):
            parFileName = "{0}.run{1}.chain{2}.par".format(fileBaseName, run, chain)
            jgsFileHandle.write("parameters in \"{0}\", chain({1})\n".format(parFileName, chain))
        jgsFileHandle.write("initialize\n")
        jgsFileHandle.write("update {0}, by({1})\n".format(warmup, thin))
        for p in pars: jgsFileHandle.write("monitor {0}, thin({1})\n".format(p, thin))
        jgsFileHandle.write("update {0}, by({1})\n".format(iter, thin))
        jgsFileHandle.write("coda *, stem({0}.run{1}.)\n".format(fileBaseName, run))
        jgsFileHandle.write("exit\n")
    return jgsFileName

def mkJagsBugFile(path, fileBaseName, model, override):
    bugFileName = os.path.join(path, "{0}.bug".format(fileBaseName))
    if not override and os.path.isfile(bugFileName):
        return bugFileName
    ## else, create the file
    with open(bugFileName, 'w') as bugFileHandle:
        bugFileHandle.write(model)
    return bugFileName

def mkRandomModelName(path):
    rid = "".join(rand.choice([s for s in string.ascii_letters], size=6))
    modelName = "anon_model_" + rid
    fileNameWildcard = os.path.join(path, "{0}.*".format(modelName))
    existingFileNames = glob.glob(fileNameWildcard)
    if len(existingFileNames) > 0:
        warnings.warn("random file name miraculously exists... trying another name")
        return mkRandomModelName(path)
    else: return modelName

## function that runs JAGS
def runJags(path, fileBaseName, run, verbose, dry_run, exp_stars, queue):
    logfilename = os.path.join(path, "{0}.run{1}.log".format(fileBaseName, run))
    logfile = open(logfilename, 'w')
    logstr = b'' ## bytes
    errfilename = os.path.join(path, "{0}.run{1}.err".format(fileBaseName, run))
    errfile = open(errfilename, 'w')
    jagsCmd = "jags"
    jgsFileName = os.path.join(path, "{0}.run{1}.jgs".format(fileBaseName, run))
    ## make a progress bar (or not)
    pbar = tqdm.tqdm(total=exp_stars) if verbose else None
    ## start JAGS
    if not dry_run:
        if verbose or queue is not None:
            proc = subprocess.Popen([jagsCmd, jgsFileName], stdout=subprocess.PIPE, stderr=errfile)
            while proc.poll() is None:
                #time.sleep(1e-3) ## TODO: needed?
                c = proc.stdout.read(1)
                logstr += c
                if c == b'*':
                    if pbar is not None:
                        pbar.update(1)
                if queue is not None:
                    queue.put(c)
            logfile.write(str(logstr, "utf-8"))
        else: ## no printing to screen, or signalling to parent process
            subprocess.call([jagsCmd, jgsFileName], stdout=logfile, stderr=errfile)
    ## write the output of JAGS to a file
    logfile.close()
    errfile.close()
    if pbar is not None:
        pbar.close()
    if queue is not None:
        queue.put(None) ## signal to monitor that job is finished
    return (logfilename, errfilename)


def containsOnlyInts(x):
    if isinstance(x, str):
        try:
            int(x)
        except ValueError:
            return False
        return True
    elif isinstance(x, list):
        return all([containsOnlyInts(e) for e in x])
    else:
        raise Exception("containsOnlyInts expects str or list")


def convertToNumeric(x, t=None):
    if t is None: ## guess the type
        if containsOnlyInts(x):
            return convertToNumeric(x, t=int)
        else:
            return convertToNumeric(x, t=float)
    else:
        if isinstance(x, str):
            return t(x)
        elif isinstance(x, list):
            return [convertToNumeric(e, t=t) for e in x]
        else:
            raise Exception("convertToNumeric expects str or list")


## function that reads and parses JAGS output
def parseJagsOutput(path, fileBaseName):
    """
    Parse the CODA output that JAGS generates.

    The result is a list of dictionaries, one for each independent chain.
    The keys of the dictionary are the parameter names.
    Vector-valued parameters are returned as follows:
    (parameter name, state index, vector index)

    Args:
        path: the directory containing the JAGS output.
        fileBaseName: an identifier used for naming files.

    Returns:
        list of dict
    """

    name_re = re.compile("[a-zA-Z0-9_\.]*")
    idx_re = re.compile("(?<=\[)[0-9,]*(?=\])")
    ## aux function for parsing parameter names and their index (for vectors)
    def getParNameAndIndex(pn):
        rn = name_re.search(pn)
        name = None if rn is None else rn.group()
        ri = idx_re.search(pn)
        idx = None if ri is None else [int(i) for i in ri.group().split(',')]
        return (name, idx)
    def getKey(pn, idx):
        return f"{pn}[{','.join([str(i) for i in idx])}]"
    ## aux function for collecting all elements of a parameter vector
    def sortTraces(pn, rawchain):
        indices = [getParNameAndIndex(key) for key in rawchain.keys()]
        indices = sorted([snd(x) for x in indices if fst(x) == pn])
        if indices[0] is None: ## scalar
            trace = rawchain[pn]
            return convertToNumeric(trace)
        elif len(indices[0]) == 1: ## 1-d array
            keys = [getKey(pn, idx) for idx in indices]
            traces = [rawchain[k] for k in keys]
            return transpose(convertToNumeric(traces))
        elif len(indices[0]) == 2: ## 2-d array
            n = np.max([idx[0] for idx in indices])
            m = np.max([idx[1] for idx in indices])
            traces = [[rawchain[getKey(pn, [i+1,j+1])] for j in range(m)] for i in range(n)]
            return np.array(convertToNumeric(traces)).transpose((2,0,1)).tolist()
        else:
            ## @todo: implement general method for n-d arrays!!
            raise Exception("n-dimensional arrays not implemented")

    ## actual function starts here...
    sams = [] ## return value (one dict for each chain)
    ## index file
    idxFileNameWildcard = os.path.join(path, "{0}.run*.index.txt".format(fileBaseName))
    idxFileNames = glob.glob(idxFileNameWildcard)
    runRegex = re.compile(os.path.join(path, "{0}.run([0-9]*).index.txt".format(fileBaseName)))
    run_ids = [runRegex.match(fn).group(1) for fn in idxFileNames]
    ## check if there are any index files
    if len(idxFileNames) == 0:
        return sams
    ## else...
    for run, idxFileName in zip(run_ids, idxFileNames):
        with open(idxFileName, 'r') as idxFileHandle:
            idxs = idxFileHandle.read().split('\n')
        idxs = [row.split(' ') for row in idxs if row != '']
        idxs = dict((row[0], (int(row[1]), int(row[2]))) for row in idxs if len(row) == 3)
        ## chain files
        chainFileNameWildcard = os.path.join(path, "{0}.run{1}.chain*.txt".format(fileBaseName, run))
        chainFileNames = glob.glob(chainFileNameWildcard)
        for chainFileName in chainFileNames:
            with open(chainFileName, 'r') as chainFileHandle:
                rawchain = chainFileHandle.read().split('\n')
            rawchain = [row.split() for row in rawchain if row != '']
            rawchain = [row[1] for row in rawchain if len(row) == 2]
            rawchain = dict((k, [rawchain[i] for i in range(fst(idxs[k])-1, snd(idxs[k]))])
                            for k in list(idxs.keys()))
            ## monitored parameter names
            monParNames = set([fst(getParNameAndIndex(key)) for key in rawchain.keys()])
            sams.append(dict((pn, sortTraces(pn, rawchain)) for pn in monParNames))
    return sams

def mergeChains(chains):
    """
    Takes a list of dictionaries, and merges them

    Args:
        chains -- a list of dicts

    Returns:
        a single dict
    """
    pars = list(chains[0].keys())
    sams = dict((p, []) for p in pars)
    for chain in chains:
        for p in pars:
            sams[p] += chain[p]
    return sams

def transposeChain(chain):
    """
    Takes a dictionary of traces and returns a list of state-dictionaries
    """
    parnames = sorted(chain.keys())
    if len(parnames) == 0:
        return []
    ## else: continue...
    n = len(chain[parnames[0]])
    states = [dict(zip(parnames, [chain[pn][i] for pn in parnames])) for i in range(n)]
    return states

## a jags model object

class JagsModel(object):
    """
    An object representing a JAGS model
    """
    def __init__(self, model_code="", file_name=None, model_name=None, path="/tmp/",
                 use_cache=False):
        self.path = path
        ## find the model code
        if file_name is None:
            self.model_code = model_code
        else:
            fileHandle = open(file_name, 'r')
            self.model_code = fileHandle.read()
            fileHandle.close()
        if model_name is None:
            self.model_name = mkRandomModelName(self.path)
            ## don't go looking for "anon_model" in e.g. the "/tmp/" folder
        else:
            self.model_name = model_name
        ## try loading the samples from disk (empy list if does not exist yet)
        if use_cache: ## try to retrieve existing chains
            self.sams = parseJagsOutput(self.path, self.model_name)
        else:
            self.sams = []
    def __str__(self): return "JAGS model '{0}'\n{1}".format(self.model_name, self.model_code)
    def __repr__(self): return self.__str__()
    def sampling(self, data={}, pars=[], chains=4, iter=1000, warmup=1000, thin=1,
                 verbose=False, dry_run=False, parallel=True):
        ## prepare files
        override = not dry_run
        ## create the directory if it does not exist
        try:
            os.makedirs(self.path)
        except OSError:
            if not os.path.isdir(self.path): raise
        jags_chains = 1 ## TODO: JAGS has an option for multiple chains too.
        jags_chain_ids = range(1, jags_chains+1) ## chain1, chain2, chain3,...
        run_ids = range(1, chains+1) ## run1, run2, run3,...
        datFileName = mkJagsDataFile(self.path, self.model_name, data, override)
        parFileNames = [mkJagsParameterFile(self.path, self.model_name, {},
                                            jags_chain, run, override)
                        for run in run_ids for jags_chain in jags_chain_ids]
        jgsFileNames = [mkJagsScriptFile(self.path, self.model_name, pars, iter, warmup,
                                         thin, jags_chains, run, override)
                        for run in run_ids]
        bugFileName = mkJagsBugFile(self.path, self.model_name, self.model_code, override)
        ## ignore 'parallel' when there is only one chain
        parallel = False if chains == 1 else parallel
        ## JAGS outputs stars... compute the number of stars for progress...
        exp_stars = (iter + warmup) // thin ## integer division
        ## make pipe for each run
        #pipes = [multiprocessing.Pipe() for _ in run_ids]
        #local_conns = [p[0] for p in pipes]
        #remote_conns = [p[1] if parallel and verbose else None for p in pipes]
        manager = multiprocessing.Manager()
        queue = manager.Queue()
        ## arguments for runJags to run the JAGS model.
        tasks = [(self.path,
                  self.model_name,
                  run,
                  verbose if not parallel else False,
                  dry_run,
                  exp_stars,
                  queue) for run in run_ids]
        if parallel:
            max_workers = max(1, multiprocessing.cpu_count()-1)
            if verbose:
                with ProcessPoolExecutor(max_workers=max_workers) as pool:
                    futures = [pool.submit(runJags, *task) for task in tasks]
                    pbar = tqdm.tqdm(total=chains*exp_stars)
                    running = chains
                    while running > 0:
                        c = queue.get()
                        if c == b'*':
                            pbar.update(1)
                        elif c is None:
                            running -= 1
                    results = [future.result() for future in futures]
            else:
                with multiprocessing.Pool(processes=max_workers) as pool:
                    results = pool.starmap(runJags, tasks)
        else:
            results = [runJags(*task) for task in tasks]
        ## show log and error files
        if verbose:
            for result in results:
                logFileHandle = open(result[0], 'r')
                errFileHandle = open(result[1], 'r')
                print("log:\n------------\n{}".format(logFileHandle.read()))
                print("errors:\n---------\n{}".format(errFileHandle.read()))
        ## parse the model output
        self.sams = parseJagsOutput(self.path, self.model_name)
        return self.sams

def calcDeltaWAIC(jm1, jm2, key1="log_like", key2=None, c1=0, c2=None, verbose=False, mass=0.95):
    """
    returns DeltaWAIC, DeltaWAIC_se, Delta_WAIC_lo, DeltaWAIC_hi
    """
    ## WAIC vec for model 1
    loglikes1 = [list(map(float, lls)) for lls in jm1.sams[c1][key1]]
    lpd_hat_vec1 = [np.log(np.mean(list(map(np.exp, lls)))) for lls in loglikes1]
    p_waic_hat_vec1 = [np.var(lls, ddof=1) for lls in loglikes1]
    elpd_waic_hat_vec1 = [x-y for x, y in zip(lpd_hat_vec1, p_waic_hat_vec1)]
    WAIC_vec1 = np.array([-2*x for x in elpd_waic_hat_vec1])
    ## WAIC vec for model 2
    if key2 is None: key2 = key1 ## assume that keys are equal
    if c2 is None: c2 = c1 ## idem for chain index
    loglikes2 = [list(map(float, lls)) for lls in jm2.sams[c2][key2]]
    lpd_hat_vec2 = [np.log(np.mean(list(map(np.exp, lls)))) for lls in loglikes2]
    p_waic_hat_vec2 = [np.var(lls, ddof=1) for lls in loglikes2]
    elpd_waic_hat_vec2 = [x-y for x, y in zip(lpd_hat_vec2, p_waic_hat_vec2)]
    WAIC_vec2 = np.array([-2*x for x in elpd_waic_hat_vec2])
    ## compute the differences
    DeltaWAIC_vec = WAIC_vec1 - WAIC_vec2
    DeltaWAIC = np.sum(DeltaWAIC_vec)
    ## number of observations
    N = len(DeltaWAIC_vec)
    DeltaWAIC_se = np.sqrt(N) * np.std(DeltaWAIC_vec)
    DeltaWAIC_lo, DeltaWAIC_hi = CImeanBootstrap(N * DeltaWAIC_vec, mass=mass)
    if verbose:
        s = "Delta WAIC = {0} (se = {1}, {2}%CI = [{3}, {4}])"
        v = (DeltaWAIC, DeltaWAIC_se, mass*100, DeltaWAIC_lo, DeltaWAIC_hi)
        print(s.format(*v))
    return (DeltaWAIC, DeltaWAIC_se, DeltaWAIC_lo, DeltaWAIC_hi)

## some plotting tools (move to another file within the module)

def density_plot(ax, xs, facecolor='purple', color='purple', alpha=1.0):
    ## TODO: sharp boundaries
    k = sts.gaussian_kde(xs)
    m = k.dataset.min()
    M = k.dataset.max()
    x = np.arange(m, M, (M-m)/100)
    v = k.evaluate(x)
    ax.fill_between(x,v,facecolor=facecolor, color=color, alpha=alpha)
    ax.plot(x, v, color=color)


def violin_plot(ax, datas, pos=None, color="violet", facecolor=None, alpha=1, midline=True):
    warnings.warn("violin plot is available from matplotlib.pyplot", DeprecationWarning)
    if facecolor is None: facecolor = color
    if pos is None: pos = list(range(len(datas)))
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for p, d in zip(pos, datas):
        if len(d) > 1:
            k = sts.gaussian_kde(d) #calculates the kernel density
            m = k.dataset.min() #lower bound of violin
            M = k.dataset.max() #upper bound of violin
            x = np.arange(m, M, (M-m)/100) # support for violin
            v = k.evaluate(x) #violin profile (density curve)
            v = v/v.max()*w #scaling the violin to the available space
            if midline:
                ax.fill_betweenx(x, p, v+p, facecolor=facecolor, color=color, alpha=alpha)
                ax.fill_betweenx(x, p, -v+p, facecolor=facecolor, color=color, alpha=alpha)
            else:
                ax.fill_betweenx(x, -v+p, v+p, facecolor=facecolor, color=color, alpha=alpha)
        elif len(d) == 1:
            ax.plot([-v+p, v+p], [d[0], d[0]], color=color, alpha=alpha)
        else:
            print("warning (violin_plot): no data in column " + str(p))
    ax.set_xlim(min(pos)-1, max(pos)+1)

## some examples/demonstrations

def papyjags_example():
    """
    give a basic example.
    """
    model = """
    model {
        for ( i in 1:N ) {
            X[i] ~ dnorm(mu,tau)
        }
        mu ~ dnorm(0,0.01)
        tau ~ dgamma(0.01,0.01)
    }
    """
    print("example...", model)
    N = 100
    data = {"N" : N, "X" : sts.norm.rvs(loc=0, scale=1, size=N)}
    pars = ["mu", "tau"]
    jm = JagsModel(model_code=model)
    print(jm)
    jm.sampling(iter=1000, warmup=1000, data=data, pars=pars, verbose=True, chains=1)
    chain = mergeChains(jm.sams)
    mus = chain["mu"]
    taus = chain["tau"]
    print("posterior means:", np.mean(mus), np.mean(taus))
    fig, (ax1, ax2) = plt.subplots(1, 2)
    density_plot(ax1, mus)
    ax1.set_xlabel("$\\mu$")
    density_plot(ax2, taus)
    ax2.set_xlabel("$\\tau$")
    fig.show()
