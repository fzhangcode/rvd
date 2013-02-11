#!/usr/bin/env python

"""rvd23.py: Compute MAP estimates for RVD2.3 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.stats import gamma
from scipy.special import gammaln, psi, polygamma

import logging


def generate_sample(phi, n=1000, nrep=1, npos=1, ncat=4):
    """Return n samples from the RVD2.3 model with nrep replciates, npos locations and ncat multinomial categories."""
    
    # Draw M Dirichlet parameters
    alpha = gamma.rvs(phi['a'], scale=phi['b'], size=ncat*npos)
    alpha = np.reshape(alpha, (npos, ncat))
    
    # Draw the sample/replicate probabilities
    theta = np.zeros((nrep, npos, ncat))
    r = np.zeros((nrep, npos, ncat), dtype=np.int64)
    for i in xrange(0, nrep):
        for j in xrange(0, npos):
            theta[i,j,:] = np.random.mtrand.dirichlet(alpha[i,:])
            r[i,j,:] = np.random.multinomial(n, theta[i,j,:])
    
    return (r, alpha, theta)

def generate_sample2(theta, n=1000, nrep=1, npos=1):
    """ Return n samples with probabilities given by theta
    """
    ncat = np.shape(theta)[0]
    
    r = np.zeros((nrep, npos, ncat), dtype=np.int64)
    for i in xrange(0, nrep):
        for j in xrange(0, npos):
            r[i,j,:] = np.random.multinomial(n, theta)

    return r
def log_cond_alpha(phi, alpha, theta):
    """Return the log conditional distribution of alpha given its 
    Markov blanket.
    alpha = ncat-dimensional vector
    theta = NxK matrix with N replicates and K categories
    """
    
    nrep, ncat =  np.shape(theta)
    logPtheta = 0;
    for i in xrange(0,nrep):
        logPtheta = logPtheta \
                    + gammaln(np.sum(alpha)) \
                    - np.sum(gammaln(alpha)) \
                    + np.sum( (alpha-1) * np.log(theta[i,:]+np.finfo(np.float).eps) )

    logPalpha = - gammaln(phi['a']) \
                - phi['a']*np.log(phi['b']) \
                + (phi['a']-1) * np.log(alpha) \
                - alpha/phi['b']
    logPalpha = np.sum(logPalpha)
    
    return (logPtheta + logPalpha)

def log_cond_theta(phi, alpha, theta, r):
    """Return the log conditional distribution of theta given its 
    Markov blanket.
    alpha = K-dimensional vector
    theta = K-dimensional vector
    r = K-dimensional vector of multinomial counts
    """
    
    n = np.sum(r) # total counts for multinomial
    
    logPr = gammaln(n+1) - np.sum(gammaln(r+1)) + np.sum(r * np.log(theta))
    
    logPtheta = gammaln(np.sum(alpha)) \
                - np.sum(gammaln(alpha)) \
                + np.sum( (alpha-1) * np.log(theta) )
                
    return (logPr + logPtheta)


def sample_alpha_mh(alpha, theta, phi, nsample=1, scale=0.5):
    """Return a sample from the posterior marginal distribution for alpha"""
    
    ncat = np.shape(alpha)[0]
    c = 0
    
    for s in xrange(0, nsample):
        # Compute the posterior likelihood for the current alpha
        logPalpha = log_cond_alpha(phi, alpha, theta)

        # Sample a non-negative value for alpha_p
        alpha_p = np.copy(alpha)
        for k in xrange(0, ncat):
            while True: 
                alpha_p[k] = np.random.normal(loc=alpha[k], scale=scale)
                if alpha_p[k] > 0: break
            
        # Compute the posterior likelihood for the proposal alpha        
        logPalpha_p = log_cond_alpha(phi, alpha_p, theta)
    
        # Compute the acceptance probability
        prA = np.exp(np.min([0.0, logPalpha_p - logPalpha]))
  
        # Return the new value of alpha
        if (np.random.rand(1) < prA): 
            alpha = np.copy(alpha_p)
            c = c + 1

    if (c/nsample) < 0.2: 
        scale = scale/2
    elif (c/nsample) > 0.8: 
        scale = scale+1
    
    return (alpha, scale)

def sample_theta_mh(alpha, theta, r, phi, nsample=1):
    """Generate a Metropilis-Hastings sample from the 
    marginal posterior of theta
    """
    
    ncat = np.shape(alpha)[0]
    
    for s in xrange(0, nsample):
        # Compute the posterior likelihood for the current theta
        logPtheta = log_cond_theta(phi, alpha, theta, r)
            
        # Sample a non-negative value for theta_p
        theta_p = np.copy(theta)
        for k in xrange(0, ncat):
            while True: 
                theta_p[k] = np.random.normal(loc=theta[k], scale=0.01)
                if 0 <= theta_p[k] <= 1: break
        theta_p = theta_p/np.sum(theta_p)
                
        # Compute the posterior likelihood for the proposal theta        
        logPtheta_p = log_cond_theta(phi, alpha, theta_p, r)
        
        # Compute the acceptance probability
        prA = np.exp(np.amin([0.0, logPtheta_p - logPtheta]))
            
        # Write the new value of alpha
        if (np.random.rand(1) < prA): 
            theta = np.copy(theta_p)
        
    return theta

def sample_theta_gibbs(alpha, r):
    """ Return a sample from the conditional, p(theta|r, alpha),
    """
    
    return np.random.mtrand.dirichlet(alpha+r)
    
def metro_gibbs(r, phi, nsample=10000):
    """Return samples from the posterior p(theta,alpha) 
    using Metropolis-within-Gibbs sampling.
    """
    
    # Get the size of the data set from r
    nrep, npos, ncat = np.shape(r)
    
    # Initialize alpha and theta
    alpha = gamma.rvs(1, scale=1, size=ncat*npos)
    alpha = np.reshape(alpha, (npos, ncat))
    
    theta = np.zeros((nrep, npos, ncat))
    for i in xrange(0, nrep):
        for j in xrange(0, npos):
            theta[i,j,:] = np.random.mtrand.dirichlet(alpha[i,:])

    # Draw nsample samples from the posterior distribution using M-H
    theta_s = np.zeros((nrep, npos, ncat, nsample))
    alpha_s = np.zeros((npos, ncat, nsample))
    sc = 0.5
    for s in xrange(0,nsample):
        for j in xrange(0,npos):
            alpha[j,:], sc = sample_alpha_mh(alpha[j,:], theta[:,j,:], phi, scale=sc)
            for i in xrange(0, nrep):
                theta[i,j,:] = sample_theta_gibbs(alpha[j,:], r[i,j,:])
                # theta[i,j,:] = sample_theta_mh(alpha[j,:], theta[i,j,:], r[i,j,:], phi)
        # Store the sample and update Gamma parameters by MLE
        alpha_s[:,:,s] = np.copy(alpha)
        theta_s[:,:,:,s] = np.copy(theta)
        # phi['a'], phi['b'] = gamma_mle(alpha_s[:,:,s])
        
    return (alpha_s, theta_s)

def gamma_mle(x):
    """ Return the maximum-likelihood estimates for shape and scale of Gamma RV.
    This uses the method proposed by T. Minka in 
    http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
    """
    
    # Initialize a
    a = 0.5 / ( np.log(np.mean(x)) - np.mean(np.log(x)) )
    
    
    # Generalized Newton updates for "a"
    while True:
        a_new = 1 / (1/a + (np.mean(np.log(x)) \
                        - np.log(np.mean(x)) \
                        + np.log(a) \
                        - psi(a) ) \
                        / 
                        (np.power(a,2) \
                         * (1/a - polygamma(1, a))))
        if (abs(a-a_new)/abs(a)) < 1e-6: break
        else: a = a_new
    
    b = np.mean(x)/a
    
    return (a, b)

if __name__ == '__main__':
    npos = 100
    ncat = 2
    nrep = 3
    nsample = 2000
    burnin = 0.2
    
    
    phi = {'a':1, 'b':20}
    # r, alpha, theta = generate_sample(phi, n=10000, npos=npos, nrep=nrep, ncat=ncat)
    
    theta = np.array([0.9, 0.1])
    ncat = np.shape(theta)
    r = generate_sample2(theta, n=1000, nrep=nrep, npos=npos)
    
    for i in xrange(0, 10):
        alpha_s, theta_s = metro_gibbs(r, phi, nsample=nsample)
    
        # Remove burn-in samples
        nps = np.int(burnin*nsample) if burnin < 1 else burnin
        nsample = nsample - nps
        alpha_s = np.delete(alpha_s, range(0, nps), 2)
        theta_s = np.delete(theta_s, nps, 3)
    
        # Estimate Gamma hyperparameters
        alpha_hat = alpha_s.mean(2)
        theta_hat = theta_s.mean(3)
        
        # Set a and b based on the reference base only
        phi['a'], phi['b'] = gamma_mle(alpha_hat[:,0])
        
        # plt.plot(np.arange(nsample), alpha_s[0,0,:], color='b')
        # plt.show()
        
        print "Estimated phi"
        print phi
    


    col = ('b', 'g', 'r', 'c')
    ind = np.arange(npos)
    
    

     
    g = alpha_hat[:,0] / np.sum(alpha_hat,1)
    plt.plot(ind, g, color=col[0], marker='x')
    for i in xrange(0,nrep):
        # plt.plot(np.arange(npos), theta[i,:,0], marker='o', ls='', color=col[i])
        plt.plot(np.arange(npos), theta_hat[i,:,0], marker='+', ls='', color=col[i])
    
    plt.show()
    