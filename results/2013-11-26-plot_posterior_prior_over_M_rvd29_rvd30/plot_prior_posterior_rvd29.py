#!/usr/bin/env python
# Plot the probability of Jeffreys' prior and the posterior Mj over a range of M.

from scipy.special import gammaln, polygamma
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import sys
import os
import h5py
import logging
import rvd29
import itertools
import math
# read depth rate = 1/100
h5Filename = './hdf5_rvd29_gibbs=20000_dc=100/Case0_1.hdf5'
loc = 104

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    M = range(1,10000)

    # Probability of Jeffreys' prior from rvd29
    (Phi, M_s, Theta, Mu, Loc, R, N) = rvd29.load_model(h5Filename)

    # mu from gibbs
    mu=np.mean(Mu, axis=0) 
    mu=np.mean(mu)
    print mu
    
    p = np.sqrt(polygamma(1, np.dot(mu,M))*(mu**2) + polygamma(1, np.dot(1-mu,M))*(1-mu)**2 -polygamma(1,M))

    M1 = M_s[loc]
    print M1

    # Use Gaussian Kernel Density Estimate for posterior Mj
    kernel = gaussian_kde(M1)  
    fig=plt.figure(figsize=(20,12))

    plt.semilogy(M,p)
    plt.semilogy(M,kernel(M).T)

    #figname = "pdf_%s" %loc
    plt.savefig('pdf.pdf')

if __name__ == '__main__':
    main()
    
