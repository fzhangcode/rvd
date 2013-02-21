#!/usr/bin/env python

"""rvd25.py: Compute MAP estimates for RVD2.5 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import logging

def main():
    n = 1000
    phi = {'mu0':0.25, 'M0':2e3, 'M':10.0}
    r, theta, mu = generate_sample(phi, n=n, N=3,  J=10)
    
    loglik = complete_ll(phi, r, n, theta, mu)
    
    phi, mu, theta = estimate_mom(r, n)
    phi, theta, mu = gibbs_map(r, n)
    plt.plot(range(0,10), mu)
    plt.plot(range(0,10), np.transpose(theta), linestyle='None', marker='+', markersize=10)
    plt.plot(range(0,10), np.tile(phi['mu0'], 10), linestyle='--', color='r')
    plt.show()

    

def generate_sample(phi, n=100, N=3, J=100):
    """Returns a from the RVD2.5 model with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    
    # Draw J location-specific error rates from a Log-Normal
    alpha0 = phi['M0']*phi['mu0']
    beta0 = phi['M0']*(1-phi['mu0'])
    mu = ss.beta.rvs(alpha0, beta0, size=J)
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*phi['M']
        beta = (1-mu[j])*phi['M']
        theta[:,j] = ss.beta.rvs(alpha, beta, size=N)
        r[:,j] = ss.binom.rvs(n, theta[:,j])
    return r, theta, mu

def complete_ll(phi, r, n, theta, mu):
    """ Return the complete data log-likelihood.
    """
    alpha0 = phi['M0']*phi['mu0'] + np.finfo(np.float).eps
    beta0 = phi['M0']*(1-phi['mu0']) + np.finfo(np.float).eps
    
    alpha = phi['M']*mu + np.finfo(np.float).eps
    beta = phi['M']*(1 - mu) + np.finfo(np.float).eps
    
    # TODO check for divide by zero in binomial logpdf when theta is close to 0 or 1
    logPmu = beta_log_pdf(mu, alpha0, beta0)
    logPtheta = beta_log_pdf(theta, alpha, beta)
    logPr = ss.binom.logpmf(r, n, theta)
    
    return np.sum(logPmu + logPtheta + logPr)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """
    
    theta = r/n
    mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps)
    
    M = np.mean( (mu*(1-mu))/(np.var(theta, 0) + 1e-6 ))
    
    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    return phi, mu, theta
    
def gibbs_map(r, n, tol=1e-4):
    """ Return MAP parameter and latent variable estimates obtained by 
    Metropolis-Hastings within Gibbs sampling.
    By default, sample 10000 M-H with a 20% burn-in. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """
    
    N, J = np.shape(r)
    
    # Initialize estimates using MoM
    phi, mu, theta = estimate_mom(r, n)
    # phi['M0'] = 2e3
    
    llCurr = complete_ll(phi, r, n, theta, mu)
    llDelta = np.inf
    # while (llDelta > tol):
    for i in xrange(0,10):
    
        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*phi['M']
        beta = (n - r) + (1-mu)*phi['M']
        theta = ss.beta.rvs(alpha+1e-3, beta+1e-3)

        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_s = sampleMuMH(theta, phi['mu0'], phi['M0'], phi['M'], mu=mu,
                          nsample=500)
        mu = np.median(mu_s, 0)
        # plt.hist(mu_s[:,0], 50, normed=1, facecolor='g', alpha=0.75)
        # plt.xlabel('mu[0]')
        # plt.ylabel('Posterior Probability')
        # plt.grid(True)
        # plt.show()
        
        # Compute MLE for parameters
        s2 = np.zeros(J)
        mu2 = (mu+np.finfo(np.float).eps)*(1-mu-np.finfo(np.float).eps)
        for j in xrange(0,J):
            s2[j] = N*np.sum(n*(theta[:,j]-mu[j])**2) / ((N-1)*N*n)
        print s2
        print mu2
        M_s = (mu2 - s2)/(s2 - (N/n)*(mu2)/N )
        phi['M'] = np.mean(M_s)
        phi['mu0'] = np.mean(mu)
        
        # Check for convergence
        llPrev = np.copy(llCurr)
        llCurr = complete_ll(phi, r, n, theta, mu)
        llDelta = (llCurr - llPrev)/np.abs(llPrev)
        print llDelta
        # llDelta = 0
    
        # Diagnosic Plots
        Brep = phi['M']/(phi['M']+n)
        print("Replicate Shrinkage Factor: %0.3f" % Brep)
        
        alpha0 = phi['mu0']*phi['M0']
        beta0 = (1-phi['mu0'])*phi['M0']
        print phi
        
    return (phi, theta, mu)

def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x)

def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), burnin=0.2, nsample=5000, thin=2):
    """ Return a sample of mu with parameters mu0 and M0.
    """
    N,J = np.shape(theta)
    alpha0 = mu0*M0 + np.finfo(np.float).eps
    beta0 = (1-mu0)*M0 + np.finfo(np.float).eps
    
    
    mu_s = np.zeros( (nsample, J) ) 
    for ns in xrange(0, nsample):
        mu_p = np.zeros(J)
        for j in xrange(0, J):
        
            # Sample from the proposal distribution
            mu_p = ss.norm.rvs(mu[j], 1e-2)
            if not (0 < mu_p < 1): continue
        
            # Log-likelihood for the proposal mu
            alpha_p = mu_p*M + np.finfo(np.float).eps
            beta_p = (1-mu_p)*M + np.finfo(np.float).eps
            logPmu_p = beta_log_pdf(mu_p, alpha0, beta0) \
                        + np.sum(beta_log_pdf(theta[:,j], alpha_p, beta_p))
                    
            # Log-likelihood for the current mu
            alpha = mu[j]*M + np.finfo(np.float).eps
            beta = (1-mu[j])*M + np.finfo(np.float).eps
            logPmu = beta_log_pdf(mu[j], alpha0, beta0) \
                        + np.sum(beta_log_pdf(theta[:,j], alpha, beta))
        
            # Accept new mu if it increases posterior pdf or by probability
            loga = logPmu_p - logPmu
            if (loga > 0 or np.log(np.random.random()) < loga): 
                mu[j] = np.copy(mu_p)
        
        # Save the new sample
        mu_s[ns,:] = np.copy(mu)
    
    mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nsample):], 0)
    mu_s = np.delete(mu_s, np.s_[::thin], 0)
    return mu_s
   
   
def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass


if __name__ == '__main__':
    main()
    