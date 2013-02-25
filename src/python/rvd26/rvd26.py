#!/usr/bin/env python

"""rvd26.py: Compute MAP estimates for RVD2.6 model."""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import logging

def main():
    n = 10
    J = 50
    phi = {'mu0':0.25, 'M0':2e3, 'a':1000, 'b':1}
    r, theta, mu, M = generate_sample(phi, n=n, J=J)
    r[:,int(J/2)] = n*np.array([0.50, 0.55, 0.45])
    loglik = complete_ll(phi, r, n, theta, mu, M)
    
    # phi, mu, theta, M = estimate_mom(r, n)
    
    phi, theta, mu, M = gibbs_map(r, n)

def generate_sample(phi, n=100, N=3, J=100):
    """Returns a sample with n reads, N replicates, and
    J locations. The parameters of the model are in the structure phi.
    """
    
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
    
    # Initialize using Stirling's approximation
    a = 0.5/(np.log(np.mean(x)) - np.mean(np.log(x)))
    b = np.mean(x)/a
    
    for i in xrange(0,5):
        inva = 1/a + (np.mean(np.log(x)) - np.log(np.mean(x)) + np.log(a) - psi(a )) \
                    / (a**2 *(1/a + polygamma(1,a)))
        a = 1/inva
        b = np.mean(x)/a
    
    return (a, b)
def gibbs_map(r, n, tol=1e-4):
    """ Return MAP parameter and latent variable estimates obtained by 
    Metropolis-Hastings within Gibbs sampling.
    By default, sample 10000 M-H with a 20% burn-in. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """
    
    if np.ndim(r) == 1: 
        J = np.shape(r)[0]
        N=1
    elif np.ndim(r) == 2:
        N, J = np.shape(r)
        
    # if n.len() == 1: n = np.tile(n, (J))
    
    # Initialize estimates using MoM
    phi, mu, theta, M = estimate_mom(r, n)
    print phi
    # phi = {'mu0':0.25, 'M0':20, 'a':1000, 'b':1}
    
    llCurr = complete_ll(phi, r, n, theta, mu, M)
    llDelta = np.inf
    # while (llDelta > tol):
    for i in xrange(0,20):
    
        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*M
        beta = (n - r) + (1-mu)*M
        theta = ss.beta.rvs(alpha+1e-3, beta+1e-3)
        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_s = sampleMuMH(theta, phi['mu0'], phi['M0'], M, mu=mu, nsample=5000)
        mu = np.median(mu_s, 0)
        muErr = np.abs(np.percentile(mu_s, (2.5, 97.5), 0) - mu)
        
        # Draw samples from p(M | a, b, theta, mu)
        M_s = sampleMMH(theta, mu, phi['a'], phi['b'], M=M, nsample=5000)
        M = np.median(M_s, 0)
        
        # Update parameter estimates
        phi['mu0'] = np.mean(mu)
        phi['M0'] = (phi['mu0']*(1-phi['mu0']))/(np.var(mu) + np.finfo(np.float).eps)
        
        phi['a'], phi['b'] = gamma_mle(M)
        print phi
        
        # Check for convergence
        llPrev = np.copy(llCurr)
        llCurr = complete_ll(phi, r, n, theta, mu, M)
        llDelta = (llCurr - llPrev)/np.abs(llPrev)
    
        # Diagnosic Plots
        fig = plt.figure()
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        
        ax1.hist(mu_s[:, int(J/2)], 50, normed=1, facecolor='g', alpha=0.75)
        ax1.set_xlabel('mu[50]')
        ax1.set_ylabel('Posterior Probability')
        ax1.xaxis.grid(True)
        ax1.yaxis.grid(True)
        
        ax2.hist(M_s[:, int(J/2)], 50, normed=1, facecolor='g', alpha=0.75)
        ax2.set_xlabel('M[50]')
        ax2.set_ylabel('Posterior Probability')
        ax1.xaxis.grid(True)
        ax1.yaxis.grid(True)

        plt.show()
        
        
        # plt.plot(range(0,J), mu, color='b')
        plt.errorbar(range(0,J), mu, yerr=muErr, color='b', linestyle='--')
        plt.plot(range(0,J), np.transpose(theta), linestyle='None', marker='+', markersize=10)
        plt.plot(range(0,J), np.transpose(r/n), linestyle='None', marker='o', markersize=6)
        plt.plot(range(0,J), np.tile(phi['mu0'], J), linestyle='--', color='r')
        plt.show()
        
    return (phi, theta, mu, M)

def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x)

def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), burnin=0.2, nsample=5000, thin=2):
    """ Return a sample of mu with parameters mu0 and M0.
    """
    if np.ndim(theta) == 1:
        N = 1
        J = np.shape(theta)[0]
    elif np.ndim(theta) > 1:
        N, J = np.shape(theta)
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
            alpha_p = mu_p*M[j] + np.finfo(np.float).eps
            beta_p = (1-mu_p)*M[j] + np.finfo(np.float).eps
            if N == 1:
                logPmu_p = beta_log_pdf(mu_p, alpha0, beta0) \
                            + np.sum(beta_log_pdf(theta[j], alpha_p, beta_p))
            elif N > 1:
                logPmu_p = beta_log_pdf(mu_p, alpha0, beta0) \
                            + np.sum(beta_log_pdf(theta[:,j], alpha_p, beta_p))
                    
            # Log-likelihood for the current mu
            alpha = mu[j]*M[j] + np.finfo(np.float).eps
            beta = (1-mu[j])*M[j] + np.finfo(np.float).eps
            if N == 1:
                logPmu = beta_log_pdf(mu[j], alpha0, beta0) \
                            + np.sum(beta_log_pdf(theta[j], alpha, beta))
            elif N > 1:
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
  
def sampleMMH(theta, mu, a, b, M=ss.gamma.rvs(1, 1), burnin=0.2, nsample=5000, thin=2):
    """ Return a sample of M with parameters a and b.
    """
    if np.ndim(theta) == 1:
        N = 1
        J = np.shape(theta)[0]
    elif np.ndim(theta) > 1:
        N, J = np.shape(theta)
       
    
    accP = np.zeros(J) # track the acceptance probability for each position
    Qsd = M/10 # keep the std, dev of the proposal
    
    M_s = np.zeros( (nsample, J) ) 
    for ns in xrange(0, nsample):
        for j in xrange(0, J):
        
            # Sample from the proposal distribution
            while True:
                M_p = ss.norm.rvs(M[j], Qsd[j])
                if M_p > 0: break;
        
            # Log-likelihood for the proposal mu
            alpha_p = mu[j]*M_p + np.finfo(np.float).eps
            beta_p = (1-mu[j])*M_p + np.finfo(np.float).eps
            logPM_p = np.sum(beta_log_pdf(mu[j], alpha_p, beta_p)) \
                        + ss.gamma.logpdf(M_p, a, scale=b)
                    
            # Log-likelihood for the current mu
            alpha = mu[j]*M[j] + np.finfo(np.float).eps
            beta = (1-mu[j])*M[j] + np.finfo(np.float).eps
            logPM = np.sum(beta_log_pdf(mu[j], alpha, beta)) \
                        + ss.gamma.logpdf(M[j], a, scale=b)

            # Accept new mu if it increases posterior pdf or by probability
            loga = logPM_p - logPM
            # if j==0: print (M_p, M[j], logPM_p, logPM)
            if (loga > 0 or np.log(np.random.random()) < loga): 
                M[j] = np.copy(M_p)
                accP[j] += 1
            
        # if ns%100 == 0:
        #     for j in xrange(0, J):
        #         if accP[j] < 0.2: Qsd[j] = 0.9*Qsd[j]
        #         elif accP[j] > 0.8: Qsd[j] = 2*Qsd[j]
        #     accP = np.zeros(J)
        
        # Save the new sample
        M_s[ns,:] = np.copy(M)
    
    M_s = np.delete(M_s, np.s_[0:np.int(burnin*nsample):], 0)
    M_s = np.delete(M_s, np.s_[::thin], 0)
    return M_s
   
def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass


if __name__ == '__main__':
    main()
    