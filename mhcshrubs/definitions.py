import os.path

## get the path to the package disrectory
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
NETMHCPAN_PATH = "/tbb/local/ii/bin" ## FIXME
## TODO JAGS path

missing_code = "missing" ## TODO: parameters
left_censored_code = "left_censored" ## the given value is an upper bound
right_censored_code = "right_censored" ## the given value is a lower bound
uncensored_code = "uncensored"

auxiliaryLowerCensBound = 0.0 ## TODO: used in JAGS. is this nescesary?

## colors used for HLA alleles depend on the locus
locusColorDict = {
    'A' : 'gray',
    'B' : 'green',
    'C' : 'purple'
}

## translate parameters in JAGS to symbols used in the manuscript
symbolDict = {
    'alpha' : "$\\mu_{\\rm spvl}$",
    'betaNodes' : "$\\alpha$",
    'betaA' : "$\\beta_A$",
    'betaB' : "$\\beta_B$",
    'betaC' : "$\\beta_C$",
    'pA' : "$p_A$",
    'pB' : "$p_B$",
    'pC' : "$p_C$",
    'tau_V' : "$\\tau_{\\rm spvl}$",
    'tau_beta' : "$\\tau_{\\rm branch}$",
    'eta' : "$\\eta$"
}
