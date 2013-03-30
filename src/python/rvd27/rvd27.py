#!/usr/bin/env python

"""rvd27.py: Compute MAP estimates for RVD2.7 model."""

import numpy as np

import scipy as sp
import scipy.stats as ss
from scipy.special import gammaln, psi, polygamma

import matplotlib.pyplot as plt
import logging
import time
import multiprocessing as mp
from itertools import repeat
import h5py
import tempfile

def main():
    import argparse
    
    # Populate our options, -h/--help is already there.
    argp = argparse.ArgumentParser(prog='rvd2', description="RVD is a hierarchical bayesian network to identify rare variants from short-read sequence data")
    argp.add_argument('--version', action='version', version='%(prog)s 2.6')
    argp.add_argument('-v', '--verbose', dest='verbose', action='count',
                    help="increase verbosity (specify multiple times for more)")
    argp.add_argument('cmd', action='store', choices=['gen', 'gibbs'])
    
    # Parse the arguments (defaults to parsing sys.argv)
    args = argp.parse_args()
    
    # TODO check what came in ton the command line and call optp.error("Useful message") to exit if all is not well

    log_level=logging.WARNING #default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >=2:
        log_level = logging.DEBUG
    
    # Set up basic configuration, out to stderr with a reasonable default format
    logging.basicConfig(level=log_level,
                        format='%(levelname)s:%(module)s:%(message)s')
                        
    # Do actual work here
    sample_run()
    
def sample_run():
    n = 1000
    J = 10
    phi = {'mu0':0.20, 'M0':2e3, 'a':1e6, 'b':1}
    r, theta, mu, M = generate_sample(phi, n=n, J=J, seedint=10)
    r[:,int(J/2)] = n*np.array([0.50, 0.55, 0.45])

    phi, theta_s, mu_s = mh_sample(r, n, 
                                        nsample=100, 
                                        thin=0, 
                                        burnin=0)
    # plot_estimate(r, n, mu_s, theta_s, phi)

def load_model(h5Filename):
    """ Returns the RVD2.6 model samples and parameters.
    Takes an hdf5 filename and returns (phi, theta_s, mu_s, M_s), where _s 
    indicates a matrix of samples where the samples are in the last dim.
    """

    h5file = h5py.File(h5Filename, 'r')
    phi = {'mu0': h5file['phi/mu0'][...],
           'M0': h5file['phi/M0'][...],
           'M': h5file['phi/M'][...]}
    theta_s = h5file['theta_s'][...]
    mu_s = h5file['mu_s'][...]
    h5file.close()
    
    return (phi, theta_s, mu_s, M)
    
def save_model(h5Filename, loc, refb, r, n, phi, theta_s, mu_s):
    """ Save the RVD2.7 model samples and parameters """
    (N, J, nsample) = np.shape(theta_s)
    
    # TODO add attributes to hdf5 file
    h5file = h5py.File(h5Filename, 'w')
    
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', data=phi['mu0'])
    h5file['phi'].create_dataset('M0', data=phi['M0'])
    h5file['phi'].create_dataset('M', data=phi['M'], chunks=True, fletcher32=True, compression='gzip')
    
    h5file.create_dataset('theta_s', data=theta_s, chunks=True, fletcher32=True, compression='gzip')
    h5file.create_dataset('mu_s', data=mu_s, chunks=True, fletcher32=True, compression='gzip')

    h5file.create_dataset('r', data=r, chunks=True, fletcher32=True, compression='gzip')
    h5file.create_dataset('n', data=n, chunks=True, fletcher32=True, compression='gzip')

    h5file.create_dataset('loc', data=loc, chunks=True, fletcher32=True, compression='gzip')
    h5file.create_dataset('refb', data=refb)

    h5file.close()
     
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
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*phi['M'][j]
        beta = (1-mu[j])*phi['M'][j]
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
    
    # Bound theta away from 0 or 1
    theta[theta < np.finfo(np.float).eps] = np.finfo(np.float).eps
    theta[theta > 1-np.finfo(np.float).eps] = 1 - np.finfo(np.float).eps
    
    logPmu = beta_log_pdf(mu, alpha0, beta0)
    logPtheta = beta_log_pdf(theta, alpha, beta)
    logPr = ss.binom.logpmf(r, n, theta)
    
    return np.sum(logPmu + logPtheta + logPr)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """
    
    theta = r/n
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps)
    
    M = (mu*(1-mu))/(np.var(theta, 0) + np.finfo(np.float).eps )

    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    return phi, mu, theta
    
def sampleLocMuMH(args):
    # Sample from the proposal distribution at a particular location
    mu, Qsd, theta, M, alpha0, beta0 = args
    
    # TODO put an escape in here to avoid an infinite loop
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
    return mu
    
def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), burnin=0, nsample=1, thin=0, pool=None):
    """ Return a sample of mu with parameters mu0 and M0.
    """
    if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
    elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
        
    alpha0 = mu0*M0 + np.finfo(np.float).eps
    beta0 = (1-mu0)*M0 + np.finfo(np.float).eps
    
    Qsd = mu/10
    mu_s = np.zeros( (nsample, J) ) 
    for ns in xrange(0, nsample):
        if pool is not None:
            args = zip(mu, Qsd, theta.T, M, repeat(alpha0, J), repeat(beta0, J))
            mu = pool.map(sampleLocMuMH, args)
        else:
            for j in xrange(0, J):
                args = (mu[j], Qsd[j], theta[:,j], M[j], alpha0, beta0)
                mu[j] = sampleLocMuMH(args)
	
        # Save the new sample
        mu_s[ns, :] = np.copy(mu)

    if burnin > 0.0:
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nsample):], 0)
    if thin > 0:
        mu_s = np.delete(mu_s, np.s_[::thin], 0)
    
    return mu_s
    
def mh_sample(r, n, nsample=10000, burnin=0.2, thin=2, pool=None):
    """ Return MAP parameter and latent variable estimates obtained by 

    Metropolis-Hastings sampling.
    By default, sample 10000 M-H with a 20% burn-in and thinning factor of 2. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """
    
    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    
    # Initialize a hdf5 file for logging model progress
    h5Filename = tempfile.NamedTemporaryFile(suffix='.hdf5').name
    logging.debug("Storing temp data in %s" % h5Filename)
    h5file = h5py.File(h5Filename, 'w')
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', (1,), dtype='f')
    h5file['phi'].create_dataset('M0', (1,), dtype='f')
    h5file['phi'].create_dataset('M', (J,), dtype='f')
    h5file.create_dataset('theta_s', (N, J, nsample), dtype='f')
    h5file.create_dataset('mu_s', (J, nsample), dtype='f')
    
    # Initialize estimates using MoM
    phi, mu, theta = estimate_mom(r, n)
    logging.debug("MoM: mu0 = %0.3e; M0 = %0.3e." % (phi['mu0'], phi['M0']) )

    # Correct MoM estimates to be non-trivial
    mu[mu < np.finfo(np.float).eps*1e4] = phi['mu0']
    theta[theta < np.finfo(np.float).eps*1e4] = phi['mu0']
    phi['M'][phi['M'] < np.finfo(np.float).eps *1e4] = 1
    
    # Sample theta, mu, M and update parameter estiamtes
    theta_s = np.zeros( (N, J, nsample) )
    mu_s = np.zeros( (J, nsample) )
    for i in xrange(0, nsample):
        if i % 100 == 0 and i > 0: logging.debug("Gibbs Iteration %d" % i)
            
        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*phi['M'] +  + np.finfo(np.float).eps
        beta = (n - r) + (1-mu)*phi['M'] + np.finfo(np.float).eps
        theta = ss.beta.rvs(alpha, beta)
        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_mh = sampleMuMH(theta, phi['mu0'], phi['M0'], phi['M'], mu=mu, nsample=50, pool=pool)
        mu = np.median(mu_mh, axis=0)
        
        # Store the sample
        theta_s[:,:,i] = np.copy(theta)
        mu_s[:,i] = np.copy(mu)
        
        # Update parameter estimates
        # phi['mu0'] = np.mean(mu)
        # phi['M0'] = (phi['mu0']*(1-phi['mu0']))/(np.var(mu) + np.finfo(np.float).eps)
        # TODO update for M
        
        # Store the current model
        h5file['phi']['mu0'][0] = phi['mu0']
        h5file['phi']['M0'][0] = phi['M0']
        h5file['phi']['M'][...] = phi['M']
        h5file['theta_s'][:,:,i] = theta
        h5file['mu_s'][:,i] = mu
        h5file.flush()
    
    # Apply the burn-in and thinning
    if burnin > 0.0:
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*nsample):], 1)
        theta_s = np.delete(theta_s, np.s_[0:np.int(burnin*nsample):], 2)
    if thin > 1:
        mu_s = np.delete(mu_s, np.s_[::thin], 1)
        theta_s = np.delete(theta_s, np.s_[::thin], 2)
    
    h5file.close()
    return (phi, theta_s, mu_s)

def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x)


def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass


if __name__ == '__main__':
    main()
    
