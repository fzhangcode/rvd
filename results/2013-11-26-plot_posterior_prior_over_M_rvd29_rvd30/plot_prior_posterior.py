from scipy.special import gammaln, polygamma
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import sys
import os
import h5py
import logging
import itertools
import math
import rvd30
import rvd29

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    locationList = (100,300,104,244)
    
    # Same experiment condition: read depth rate = 1/100
    h5Filename_rvd29 = './hdf5_rvd29_gibbs=20000_dc=100/Case10_0.hdf5'
    h5Filename_rvd30 = './hdf5_rvd30_gibbs=20000_dc=100/Case10_0.hdf5'
    M = range(1,10000)
    fig=plt.figure(figsize=(20, 12))

    """ Probability of Jeffreys' prior from rvd29"""
    logging.debug("Plot rvd29 prior and posterior")
    for loc in locationList:
        logging.debug("Processing location: %d" % loc)
        ax = fig.add_subplot(2,4,locationList.index(loc)+1)
        (Phi, M_s, Theta, Mu, Loc, R, N) = rvd29.load_model(h5Filename_rvd29)
        mu=np.mean(Mu, axis=0) 
        mu=np.mean(mu)
        p = np.sqrt(polygamma(1, np.dot(mu,M))*(mu**2) + polygamma(1, np.dot(1-mu,M))*(1-mu)**2 -polygamma(1,M))
        sum = np.sum(p)
        print (sum)
        M1 = M_s[loc]       
        # Use Gaussian Kernel Density Estimate for posterior Mj
        kernel = gaussian_kde(M1) 
        norm = kernel.integrate_box(0,10000)
        print (norm)
        ax.semilogy(M, p/sum)
        ax.semilogy(M, kernel(M).T/norm)
        ax.legend( ['prior','posterior'],loc=9, ncol=9)
        subtitle='position='+str(loc)
        ax.set_title(subtitle)
        if locationList.index(loc)==0:
           ax.set_ylabel('Probability of Jeffreys prior',fontsize=15)

    """ Probability of lognormal prior from rvd30"""
    logging.debug("Plot rvd30 prior and posterior")
    for loc in locationList:
        logging.debug("Processing location: %d" % loc)
        ax = fig.add_subplot(2,4,locationList.index(loc)+len(locationList)+1)
        (Phi, M_s, Theta, Mu, Loc, R, N) = rvd30.load_model(h5Filename_rvd30)
        mu=np.mean(Mu, axis=0) 
        mu=np.mean(mu)
        p = ss.lognorm.pdf(M, Phi['delta'], scale=np.exp(Phi['mu_M']))
        sum = np.sum(p)
        print (sum)
        M1 = M_s[loc]
        # Use Gaussian Kernel Density Estimate for posterior Mj
        kernel = gaussian_kde(M1)
        norm = kernel.integrate_box(0,10000)		
        ax.semilogy(M, p/sum)
        ax.semilogy(M, kernel(M).T/norm)
        ax.set_title('location %d' % loc, fontsize=15)
        subtitle='position='+str(loc)
        ax.set_title(subtitle)
        ax.legend( ['prior','posterior'],loc=8, ncol=8)
        subtitle='position='+str(loc)
        ax.set_title(subtitle)
        ax.set_xlabel('M value',fontsize=15)
        if locationList.index(loc)==0:
           ax.set_ylabel('Probability of log-normal prior',fontsize=15)
        #ax.set_xlim((-0.03,1.03))
        #ax.set_ylim((-0.03,1.03))

    title='post_prior_10'
    figformat='.pdf'
    plt.savefig(title+figformat)

if __name__ == '__main__':
    main()

