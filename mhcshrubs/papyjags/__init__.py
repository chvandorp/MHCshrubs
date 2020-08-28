from .papyjags import * ## TODO: export only what is needed!
import warnings

if not isJagsInstalled():
    warnings.warn("JAGS is not installed or not in PATH")
