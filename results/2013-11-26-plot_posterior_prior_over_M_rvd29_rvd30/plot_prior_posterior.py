from scipy.special import gammaln, polygamma
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
import sys
import os
import h5py
import logging
import math
import rvd30
import rvd29

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    # Actual positions are 200,300,205,265
    locationList = (199,299,184,264)
    
    # Same experiment condition: read depth rate = 1/100
    h5Filename_rvd29 = 'Case10_jeffreys.hdf5'
    h5Filename_rvd30 = 'Case10_lognormal.hdf5'
    
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
        M1 = M_s[loc]       
    # Use Gaussian Kernel Density Estimate for posterior Mj
        kernel = gaussian_kde(M1) 
        norm = kernel.integrate_box(0,10000)
        ax.semilogy(M, p/sum, 'r--')
        ax.semilogy(M, kernel(M).T/norm) 
        ax.legend( ['prior','posterior'],loc=1)
        subtitle='position='+str(loc+1)
        ax.set_title(subtitle)
        if locationList.index(loc)==0:
           ax.set_ylabel('Jeffrey\'s prior',fontsize=15)
        ax.set_ylim(0.000001,0.01)
        ax.yaxis.grid()
        ax.xaxis.grid()

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
        M1 = M_s[loc]
    # Use Gaussian Kernel Density Estimate for posterior Mj
        kernel = gaussian_kde(M1)
        norm = kernel.integrate_box(0,10000)		
        ax.semilogy(M, p/sum,'r--')
        ax.semilogy(M, kernel(M).T/norm) 
        ax.set_title('location %d' % loc, fontsize=15)
        subtitle='position='+str(loc)
        ax.set_title(subtitle)
        ax.legend( ['prior','posterior'],loc=1)
        subtitle='position='+str(loc+1)
        ax.set_title(subtitle)
        ax.set_xlabel('M value',fontsize=15)
        if locationList.index(loc)==0:
           ax.set_ylabel('log-normal prior',fontsize=15)
        ax.set_ylim(0.000001,0.01)
        ax.yaxis.grid()
        ax.xaxis.grid()
    plt.show()
    #plt.savefig('post_prior_10.pdf')

if __name__ == '__main__':
    main()

