import numpy as np
from scipy.special import erf as erf

def normal_prior(par,p):
    mu = p[0]
    sigma = p[1]
    return -0.5*(par-mu)**2/sigma**2


def uniform_prior(par,p):
    lo = p[0]
    hi = p[1]
    if lo <= par <= hi:
        return 0.0
    return -np.inf


def truncated_normal_prior(par,p):
    """
    mu = p[0]
    sigma = p[1]
    lo = p[2]
    hi = p[3]
    """
    mu = p[0]
    sigma = p[1]
    lo = p[2]
    hi = p[3]
    if lo <= par <= hi:
        return -0.5*(par-mu)**2/sigma**2 - np.log(sigma) - np.log(erf((hi-mu)/(np.sqrt(2)*sigma)) - erf((lo-mu)/(np.sqrt(2)*sigma)))
    return -np.inf
