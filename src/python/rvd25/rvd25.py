#!/usr/bin/env python

"""rvd25.py: Compute MAP estimates for RVD2.5 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import logging

def main():
    phi = {'mu0':0.1, 'M0':200.0, 'a':1000, 'b':1}
    r, theta, mu, M = generate_sample(phi, n=100, J=10)
    print r
    
    loglik = complete_ll(phi, r, 100, theta, mu, M)
    print loglik
    

def generate_sample(phi, n=100, N=3, J=100):
    """Returns a from the RVD2.5 model with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    
    # Draw J location-specific error rates from a Log-Normal
    alpha0 = phi['M0']*phi['mu0']
    beta0 = phi['M0']*(1-phi['mu0'])
    mu = ss.beta.rvs(alpha0, beta0, size=J)
    
    # Draw precision parameter M from Gamma
    M = ss.gamma.rvs(phi['a'], scale=phi['b'])
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*M
        beta = (1-mu[j])*M
        theta[:,j] = ss.beta.rvs(alpha, beta, size=N)
        r[:,j] = ss.binom.rvs(n, theta[:,j])
    return r, theta, mu, M

def complete_ll(phi, r, n, theta, mu, M):
    """ Return the complete data log-likelihood.
    """
    alpha0 = phi['M0']*phi['mu0']
    beta0 = phi['M0']*(1-phi['mu0'])
    
    alpha = mu*M
    beta = (1-mu)*M
    
    # TODO check for divide by zero in binomial logpdf when theta is close to 0 or 1
    logPmu = ss.beta.logpdf(mu, alpha0, beta0)
    logPM = ss.gamma.logpdf(M, phi['a'], phi['b'])
    logPtheta = ss.beta.logpdf(theta, alpha, beta)
    logPr = ss.binom.logpmf(r, n, theta)
    
    return np.sum(logPmu + logPM + logPtheta + logPr)

def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass


if __name__ == '__main__':
    main()
    