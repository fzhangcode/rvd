#!/usr/bin/env python
# Plot the probability of lognormal prior and the posterior Mj over a range of M.

from scipy.special import gammaln, polygamma
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import sys
import os
import h5py
import logging
import rvd30
import itertools
import math
# read depth rate = 1/100
h5Filename = '../2013-11-23_test_rvd30_synthetic_data/100/Case0_1.hdf5'
loc = 0

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    M = range(1,10000)

    # Probability of lognormal prior from rvd30
    (Phi, M_s, Theta, Mu, Loc, R, N) = rvd30.load_model(h5Filename)

    # mu from gibbs
    mu=np.mean(Mu, axis=0) 
    mu=np.mean(mu)

    print mu
    print Phi['delta']
    print Phi['mu_M']

    p = ss.lognorm.pdf(M, Phi['delta'], scale=np.exp(Phi['mu_M']))
    print p

    M1 = M_s[loc]
    
    # Use Gaussian Kernel Density Estimate for posterior Mj
    kernel = gaussian_kde(M1)  
    fig=plt.figure(figsize=(20,12))

    plt.semilogy(M,p)
    plt.semilogy(M,kernel(M).T)

    #figname = "pdf_%s" %loc
    plt.savefig('pdf.pdf')

if __name__ == '__main__':
    main()
    
