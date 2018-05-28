from __future__ import (print_function, division)
from builtins import (map, str, zip, range)
import pystan
import pickle
import hashlib
import os

## exported names

normal_prior_code = 0
laplace_prior_code = 1

## compile and cache, ore retreive a cached model

def CachedStanModel(model_code, model_name=None, path=None, **kwargs):
    """Adapted from the pystan docs"""
    code_hash = hashlib.md5(model_code.encode('ascii')).hexdigest()
    if model_name is None:
        cache_fn = 'cached-model-{}.pkl'.format(code_hash)
    else:
        cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
    if path is not None:
        cache_fn = os.path.join(path, cache_fn)
        ## if the directory path does not exist: try to make it
        try:
            os.makedirs(path)
        except OSError:
            if not os.path.isdir(path): raise
    try:
        sm = pickle.load(open(cache_fn, 'rb'))
    except:
        sm = pystan.StanModel(model_code=model_code, **kwargs)
        with open(cache_fn, 'wb') as f: pickle.dump(sm, f)
    return sm
