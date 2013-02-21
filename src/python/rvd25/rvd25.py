#!/usr/bin/env python

"""rvd25.py: Compute MAP estimates for RVD2.5 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import logging


def generate_sample(phi, n=1000, N=1, J=1):
    """Return n samples from the RVD2.5 model with N replicates, 
    J locations. The parameters of the model are in the structure phi.
    """
    
    # Draw J location-specific error rates from a Log-Normal
    mu = np.zeros(J)
    
def complete_ll(phi, r, n, theta, mu, M):
    """ Return the complete data log-likelihood.
    """
    alpha = mu*M
    beta = (1-mu)*M
    
    logPmu = lognormal(phi['mu0'], phi['sigma0'])
    logPM = gamma(phi['a'], phi['b'])
    logPtheta = beta(alpha, beta)
    logPr = binomial(r, theta, n)
    
    return logPmu + logPM + logPtheta + logPr
    

def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass
    

def main():
    pass


if __name__ == '__main__':
    main()
    