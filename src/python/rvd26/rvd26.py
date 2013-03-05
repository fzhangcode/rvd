#!/usr/bin/env python

"""rvd26.py: Compute MAP estimates for RVD2.6 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import logging

import multiprocessing as mp

def main():
    pool = mp.Pool(processes=4) # start 4 worker processes
    n = 10
    J = 5
    phi = {'mu0':0.25, 'M0':2e3, 'a':1000, 'b':1}
    r, theta, mu, M = generate_sample(phi, n=n, J=J)
    r[:,int(J/2)] = n*np.array([0.50, 0.55, 0.45])
    loglik = complete_ll(phi, r, n, theta, mu, M)
    
    phi, theta_s, mu_s, M_s = mh_sample(r, n, nsample=5)
    plot_estimate(r, n, mu_s, theta_s, phi)

def plot_estimate(r, n, mu_s, theta_s, phi):
    
    mu = np.median(mu_s, 1)
    theta = np.median(theta_s, 2)
    
    (N, J) = np.shape(theta)
    
    muErr = np.abs(np.percentile(mu_s, (2.5, 97.5), 1) - mu)
    # plt.plot(range(0,J), mu, color='b')

    plt.errorbar(np.arange(J)+1, mu, yerr=muErr, color='b', linestyle='--')
    plt.plot(np.arange(J)+1, np.transpose(theta), linestyle='None', marker='+', markersize=10)
    plt.plot(np.arange(J)+1, np.transpose(r/n), linestyle='None', marker='o', markersize=6, alpha=0.5)
    plt.plot(np.arange(J)+1, np.tile(phi['mu0'], J), linestyle='--', color='r')
    
    plt.title("Error Rate Estimate using RVD2.6")
    plt.xlim((0, J+1))
    plt.xlabel('Location')
    plt.ylabel('Error Rate Estimate')
    plt.show()

def generate_sample(phi, n=100, N=3, J=100, seedint=None):
    """Returns a sample with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    
    if seedint is not None: np.random.seed(seedint)
    
    # Draw J location-specific error rates from a Beta
    alpha0 = phi['M0']*phi['mu0']
    beta0 = phi['M0']*(1-phi['mu0'])
    mu = ss.beta.rvs(alpha0, beta0, size=J)
    
    # Draw J location-specific precisions from a Gamma
    M = ss.gamma.rvs(phi['a'], scale=phi['b'], size=J)
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*M[j]
        beta = (1-mu[j])*M[j]
        theta[:,j] = ss.beta.rvs(alpha, beta, size=N)
        r[:,j] = ss.binom.rvs(n, theta[:,j])
    return r, theta, mu, M

def complete_ll(phi, r, n, theta, mu, M):
    """ Return the complete data log-likelihood.
    """
    alpha0 = phi['M0']*phi['mu0'] + np.finfo(np.float).eps
    beta0 = phi['M0']*(1-phi['mu0']) + np.finfo(np.float).eps
    
    alpha = M*mu + np.finfo(np.float).eps
    beta = M*(1 - mu) + np.finfo(np.float).eps
    
    # Bound theta away from 0 or 1
    theta[theta < np.finfo(np.float).eps] = np.finfo(np.float).eps
    theta[theta > 1-np.finfo(np.float).eps] = 1 - np.finfo(np.float).eps
    
    logPM = ss.gamma.logpdf(M, phi['a'], scale=phi['b'])
    logPmu = beta_log_pdf(mu, alpha0, beta0)
    logPtheta = beta_log_pdf(theta, alpha, beta)
    logPr = ss.binom.logpmf(r, n, theta)
    
    return np.sum(logPmu + logPM + logPtheta + logPr)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """
    
    theta = r/n
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps)
    
    M = (mu*(1-mu))/(np.var(theta, 0) + np.finfo(np.float).eps )
    
    b = np.var(M)/np.mean(M)
    a = np.mean(M)/b
    
    phi = {'mu0':mu0, 'M0':M0, 'a':a, 'b':b}
    return phi, mu, theta, M
    
def gamma_mle(x):
    """ Return the maximum likelihood shape and scale parameters
    """
    
    if x.size == 0:
        return (np.nan, np.nan)
    
    # Initialize using Stirling's approximation
    a = 0.5/(np.log(np.mean(x)) - np.mean(np.log(x)))
    b = np.mean(x)/a
    
    for i in xrange(0,5):
        inva = 1/a + (np.mean(np.log(x)) - np.log(np.mean(x)) + np.log(a) - psi(a )) \
                    / (a**2 *(1/a + polygamma(1,a)))
        a = 1/inva
        b = np.mean(x)/a
    
    return (a, b)
def mh_sample(r, n, tol=1e-4, nsample=5000, burnin=0.2, thin=2):
    """ Return MAP parameter and latent variable estimates obtained by 
    Metropolis-Hastings sampling.
    By default, sample 10000 M-H with a 20% burn-in. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """
    
    def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), burnin=0, nsample=1, thin=0):
        """ Return a sample of mu with parameters mu0 and M0.
        """
        def sampleLocMuMH(mu, Qsd, theta, M, alpha0, beta0, out_q=None):
            # Sample from the proposal distribution at a particular location
            while True:
                mu_p = ss.norm.rvs(mu, Qsd)
                if (0 < mu_p < 1): break
        
            # Log-likelihood for the proposal mu
            alpha_p = mu_p*M + np.finfo(np.float).eps
            beta_p = (1-mu_p)*M + np.finfo(np.float).eps
            logPmu_p = beta_log_pdf(mu_p, alpha0, beta0) \
                        + np.sum(beta_log_pdf(theta, alpha_p, beta_p))
                    
            # Log-likelihood for the current mu
            alpha = mu*M + np.finfo(np.float).eps
            beta = (1-mu)*M + np.finfo(np.float).eps
            logPmu = beta_log_pdf(mu, alpha0, beta0) \
                        + np.sum(beta_log_pdf(theta, alpha, beta))
                            
            # Accept new mu if it increases posterior pdf or by probability
            loga = logPmu_p - logPmu
            if (loga > 0 or np.log(np.random.random()) < loga): 
                mu = mu_p
        
            if out_q is None: return mu
            else: out_q.put(mu)
    
        if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
        elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
        
        alpha0 = mu0*M0 + np.finfo(np.float).eps
        beta0 = (1-mu0)*M0 + np.finfo(np.float).eps
    
        Qsd = mu/10
        mu_s = np.zeros( (nsample, J) ) 
        for ns in xrange(0, nsample):
            pool = mp.Pool(processes=4)
            for j in xrange(0, 1):
                # mu[j] = sampleLocMuMH(mu[j], Qsd[j], theta[:,j], M[j], alpha0, beta0)
                result = pool.apply_async(sampleLocMuMH, (mu[j], Qsd[j], theta[:,j], M[j], alpha0, beta0))
            # Save the new sample
            mu_s[ns, :] = np.copy(mu)
    
        if burnin > 0.0:
            mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nsample):], 0)
        if thin > 0:
            mu_s = np.delete(mu_s, np.s_[::thin], 0)
    
        return mu_s



    def sampleMMH(theta, mu, a, b, M=ss.gamma.rvs(1, 1), burnin=0.0, nsample=1, thin=0):
        """ Return a sample of M with parameters a and b.
        """
        def sampleLocMMH(M, Qsd, theta, mu, a, b, out_q=None):
            # Sample from the proposal distribution
            while True:
                M_p = ss.norm.rvs(M, Qsd)
                if M_p > 0: break
        
            # Log-likelihood for the proposal mu
            alpha_p = mu*M_p + np.finfo(np.float).eps
            beta_p = (1-mu)*M_p + np.finfo(np.float).eps
            logPM_p = np.sum(beta_log_pdf(mu, alpha_p, beta_p)) \
                        + ss.gamma.logpdf(M_p, a, scale=b)
                    
            # Log-likelihood for the current mu
            alpha = mu*M + np.finfo(np.float).eps
            beta = (1-mu)*M + np.finfo(np.float).eps
            logPM = np.sum(beta_log_pdf(mu, alpha, beta)) \
                        + ss.gamma.logpdf(M, a, scale=b)

            # Accept new mu if it increases posterior pdf or by probability
            loga = logPM_p - logPM
            # if j==0: print (M_p, M[j], logPM_p, logPM)
            if (loga > 0 or np.log(np.random.random()) < loga): 
                M = np.copy(M_p)
            if out_q is None: return M
            else: out_q.put(M)
        
        if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
        elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
       
        Qsd = M/10 # keep the std, dev of the proposal
    
        M_s = np.zeros( (nsample, J) ) 
        for ns in xrange(0, nsample):
            for j in xrange(0, J):
                M[j] = sampleLocMMH(M[j], Qsd[j], theta[:,j], mu[j], a, b)
            M_s[ns,:] = np.copy(M)
    
        if burnin > 0.0:
            M_s = np.delete(M_s, np.s_[0:np.int(burnin*nsample):], 0)
        if thin > 0:
            M_s = np.delete(M_s, np.s_[::thin], 0)
        return M_s
    
    
    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    
    # Initialize estimates using MoM
    phi, mu, theta, M = estimate_mom(r, n)
    logging.debug("MoM Parameter Estimate")
    logging.debug(phi)

    theta_s = np.zeros( (N, J, nsample) )
    mu_s = np.zeros( (J, nsample) )
    M_s = np.zeros( (J, nsample) )
    for i in xrange(0, nsample):
        if i % 100 == 0: print "Sampling epoch: %d" % i

        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*M
        beta = (n - r) + (1-mu)*M
        theta = ss.beta.rvs(alpha, beta)
        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_mh = sampleMuMH(theta, phi['mu0'], phi['M0'], M, mu=mu, burnin=0.2, nsample=50)
        mu = np.median(mu_mh, 0)
        
        # Draw samples from p(M | a, b, theta, mu)
        M_mh = sampleMMH(theta, mu, phi['a'], phi['b'], M=M, burnin=0.2, nsample=50)
        M = np.median(M_mh, 0)
        
        # Update parameter estimates
        phi['mu0'] = np.mean(mu)
        phi['M0'] = (phi['mu0']*(1-phi['mu0']))/(np.var(mu) + np.finfo(np.float).eps)
        phi['a'], phi['b'] = gamma_mle(M)
    
        # Store the sample
        theta_s[:,:,i] = np.copy(theta)
        mu_s[:,i] = np.copy(mu)
        M_s[:,i] = np.copy(M)
    
    # Apply the burn-in and thinning
    M_s = np.delete(M_s, np.s_[0:np.int(burnin*nsample):], 1)
    M_s = np.delete(M_s, np.s_[::thin], 1)    
    mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nsample):], 1)
    mu_s = np.delete(mu_s, np.s_[::thin], 1)
    theta_s = np.delete(theta_s, np.s_[0:np.int(burnin*nsample):], 2)
    theta_s = np.delete(theta_s, np.s_[::thin], 2)
        
    return (phi, theta_s, mu_s, M_s)

def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x)


def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass


def oddsratio():
    pass
if __name__ == '__main__':
    main()
    
