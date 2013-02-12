#!/usr/bin/env python

"""rvd24.py: Compute MAP estimates for RVD2.3 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.stats import gamma
from scipy.special import gammaln, psi, polygamma

import logging

def main():
    npos = 10
    ncat = 2
    nrep = 4
    nsample = 50000
    burnin = 0.2
    
    cv = 0.01
    
    
    phi = {'a':np.array([np.sqrt(1/cv), np.sqrt(1/cv)]), 
           'b':np.array([1*np.sqrt(cv), 9*np.sqrt(cv)]), 
           'm':np.array([0.0, 0.0])}
    # r, alpha, theta = generate_sample(phi, n=1000, npos=npos, nrep=nrep)
    
    # Generate more realistic samples
    theta = np.array([0.1, 0.9])
    ncat = np.shape(theta)[0]
    r = generate_sample2(theta, n=100, nrep=nrep, npos=npos)
    
    for i in xrange(0, 10):
        alpha_s, theta_s = metro_gibbs(r, phi, nsample=nsample)
    
        # Remove burn-in samples
        nps = np.int(burnin*nsample) if burnin < 1 else burnin
        nsample = nsample - nps
        alpha_s = np.delete(alpha_s, range(0, nps), 2)
        theta_s = np.delete(theta_s, nps, 3)
    
        alpha_hat = alpha_s.mean(2)
        theta_hat = theta_s.mean(3)
        
        # Fit the gamma parameters
        for k in xrange(0, ncat):
            # phi['a'][k], phi['b'][k] = gamma_mle(alpha_hat[:,k])
            phi['a'][k], fit_loc, phi['b'][k] = ss.gamma.fit(alpha_hat[:,k])
        
        for k in xrange(0,ncat):
            rv = gamma(phi['a'][k], scale=phi['b'][k])
            x = np.linspace(0, np.minimum(rv.dist.b, 10), num=1000)
            plt.plot(x,rv.pdf(x))
        plt.show()
        # plt.plot(np.arange(npos), alpha_hat[:,0], color='b')
        # plt.show()
        
        
        print "Estimated phi"
        print phi
    


    col = ('b', 'g', 'r', 'c')
    ind = np.arange(npos)
     
    g_hat = alpha_hat[:,0] / np.sum(alpha_hat,1)
    plt.plot(ind, g_hat, color=col[0], marker='x')
    
    # g = alpha[:,0] / np.sum(alpha,1)
    # plt.plot(ind, g, color=col[0], marker='o', ls='-')

    for i in xrange(0,nrep):
        # plt.plot(np.arange(npos), theta[i,:,0], marker='o', ls='', color=col[i])
        plt.plot(np.arange(npos), theta_hat[i,:,0], marker='+', ls='')
    
    plt.show()

def generate_sample(phi, n=1000, nrep=1, npos=1):
    """Return n samples from the RVD2.4 model with nrep replciates, 
    npos locations and ncat multinomial categories.
    """
    
    ncat = np.shape(phi['a'])[0]
    
    # Draw M Dirichlet parameters
    alpha = np.zeros((npos, ncat))
    for k in xrange(0, ncat):
        alpha[:,k] = gamma.rvs(phi['a'][k], scale=phi['b'][k], size=npos)
    
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
def log_cond_alpha(a, b, alpha, theta):
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

                    
    logPalpha = - gammaln(a) \
                - a*np.log(b) \
                + (a-1) * np.log(alpha) \
                - alpha/b
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


def sample_alpha_mh(alpha, theta, a, b, nsample=1, scale=0.5):
    """Return a sample from the posterior marginal distribution for alpha"""
    
    ncat = np.shape(alpha)[0]
    c = 0
    
    for s in xrange(0, nsample):
        # Compute the posterior likelihood for the current alpha
        logPalpha = log_cond_alpha(a, b, alpha, theta)

        # Sample a non-negative value for alpha_p
        alpha_p = np.copy(alpha)
        for k in xrange(0, ncat):
            while True: 
                alpha_p[k] = np.random.normal(loc=alpha[k], scale=scale)
                if alpha_p[k] > 0: break
            
        # Compute the posterior likelihood for the proposal alpha        
        logPalpha_p = log_cond_alpha(a, b, alpha_p, theta)
    
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
    
    

def estimate_alpha(theta):
    """ Estimate alpha by {newton-raphson} using only theta
    Input: theta is an NxK matrix with N K-dimensional observations
    """
    N, K = np.shape(theta)

    log_theta_hat = np.mean(np.log(theta),0) # observed sufficient statistics
    
    # initialize alpha
    alpha = np.ones(K)
    
    cnt = 0
    while True:
        alpha0 = np.sum(alpha) # Dirichlet concentration (precision) parameter
 
        # log-likelihood p(theta|alpha)
        F = N * ( gammaln(alpha0) - np.sum(gammaln(alpha)) \
                  + np.sum((alpha-1)*log_theta_hat) )
    
        # derivative of the log-likelihood p(theta|alpha)
        delF = np.zeros(K)
        for k in xrange(0, K):
            delF[k] = N * (psi(alpha0) - psi(alpha[k]) + log_theta_hat[k])
    
        # hessian of the log-likelihood p(theta|alpha)
        q = -N * polygamma(1, alpha)
        c = N * polygamma(1, alpha0)
        C = np.zeros( (K,K) )
        C.fill(c)
        H = np.diagflat(q) + C
    
        # update alpha estimate by Newton-Raphson Method
        alpha_old = np.copy(alpha)
        
        la = 1.0
        while True:
            alpha = alpha_old - la * np.dot(np.linalg.pinv(H), delF)
            if np.any(alpha <= 0): la = la/2
            else: break
        cnt = cnt+1
        
        # assess convergence
        del_alpha = np.amax( np.abs(alpha-alpha_old)/np.abs(alpha_old) )
        if (del_alpha < 1e-3 or cnt > 50): break
    return alpha
    

def metro_gibbs(r, n=1000, burn=0.2, thin=2):
    """Return samples from the posterior p(theta,alpha) 
    using Metropolis-within-Gibbs sampling.
    """
    
    # Get the size of the data set from r
    N, M, K = np.shape(r)
    
    # Initialize alpha and theta using IPFP
    a_hat, b_hat, alpha, theta = ipfp(r)
    print "IPFP Done."
    print alpha

    # Draw nsample samples from the posterior distribution using M-H
    theta_s = np.zeros((N, M, K, n))
    alpha_s = np.zeros((M, K, n))
    sc = 0.5
    for s in xrange(0, n):
        for j in xrange(0, M):
            alpha[j,:], sc = sample_alpha_mh(alpha[j,:], 
                                             theta[:,j,:], 
                                             a_hat, b_hat, 
                                             scale=0.5)
            for i in xrange(0, N):
                theta[i,j,:] = sample_theta_gibbs(alpha[j,:], r[i,j,:])
                
        # Store the sample (and update Gamma parameters by MLE)
        alpha_s[:,:,s] = np.copy(alpha)
        theta_s[:,:,:,s] = np.copy(theta)
    
    alpha_hat = np.mean(alpha_s,2)
    theta_hat = np.mean(theta_s,3)
    
    a_hat, b_hat = gamma_mle(alpha_hat)
        
    return (a_hat, b_hat, alpha_hat, theta_hat)

def gamma_mle(x):
    """ Return the maximum-likelihood estimates for shape and scale of Gamma RV.
    This uses the method proposed by T. Minka in 
    http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
    """
    
    npos =  np.shape(x)[0]
    
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



def gamma_mom(x):
    """ Return the method-of-moments estimates for shape and scale of Gamma RV.
    This uses the method proposed by T. Minka in 
    http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
    """
    
    M = np.shape(x)[0]
    b = np.var(x) / np.mean(x)
    a = np.mean(x) / b
    
    return (a, b)
    
def ipfp(r, n=1000, T=10):
    
    # get the size of the data set from r
    N, M, K = np.shape(r)
    
    # initialize alpha
    alpha = np.ones( (M,K) )
    
    for t in xrange(0, T):
        # draw n samples from the posterior distribution p(theta|r, alpha)
        theta_s = np.zeros((N, M, K, n))
        for s in xrange(0, n):
            for j in xrange(0, M):
                for i in xrange(0, N):
                    theta_s[i, j, :, s] = sample_theta_gibbs( alpha[j, :], 
                                                              r[i, j, :] )
        # return point estimate from posterior p(theta|r, alpha)
        theta = np.mean(theta_s, 3)
    
        # estimate alpha | theta
        for j in xrange(0, M):
            alpha[j,:] = estimate_alpha(theta[:,j,:])
    
        # estimate a,b | alpha
        a = np.zeros(K)
        b = np.zeros(K)
        for k in xrange(0, K):
            # a[k], b[k] = gamma_mle(alpha[:,k])
            a[k], b[k] = gamma_mom(alpha[:,k])
    return (a, b, alpha, theta)

def gamma_like(p, x, sign=1.0):
    a = p[0,:]
    b = p[1,:]
    
    ll = -a * np.log(b) - np.gammaln(a) + (a - 1) * np.log(x) - x/b
    ll = sign * np.sum(ll)
    
    return ll
    

if __name__ == '__main__':
    M = 100 # number of locations
    K = 2 # number of categories
    N = 3 # number of replicates
    
    # cv = 0.01 # coefficient of variation of alpha
    # phi = {'a':np.array([np.sqrt(1/cv), np.sqrt(1/cv)]), 
    #        'b':np.array([1*np.sqrt(cv), 9*np.sqrt(cv)])}
    # r, alpha, theta = generate_sample(phi, n=1000, npos=M, nrep=N)
    # alpha0 = np.sum(alpha, 1)
    
    # Generate more realistic samples
    theta = np.array([0.1, 0.9])
    r = generate_sample2(theta, n=5, nrep=N, npos=M)
    
    a_hat, b_hat, alpha_hat, theta_hat = ipfp(r)
    alpha0_hat = np.sum(alpha_hat,1)
    
    print (a_hat, b_hat)
    
    
    # a_hat, b_hat, alpha_hat, theta_hat = metro_gibbs(r, n=1000)
    # alpha0_hat = np.sum(alpha_hat,1)
    # print a_hat
    # print b_hat
    
    g_hat = np.copy(alpha_hat)
    for k in xrange(0, K):
        g_hat[:,k] = g_hat[:,k] / alpha0_hat
        
    # g = np.copy(alpha)
    # for k in xrange(0, K):
    #     g[:,k] = g[:,k] / alpha0
    #     
    # plt.plot(np.arange(M), g, marker='^', ls='-')
    plt.plot(np.arange(M), g_hat, marker='o', ls='-')
    
    # plt.plot(np.arange(M), np.transpose(theta[:,:,0]), marker='o', ls='')
    plt.plot(np.arange(M), np.transpose(theta_hat[:,:,0]), marker='^', ls='')
    
    plt.show()
    