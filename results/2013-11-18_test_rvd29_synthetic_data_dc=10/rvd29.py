#!/usr/bin/env python


from __future__ import print_function
from __future__ import division

import numpy as np

import scipy.stats as ss
from scipy.special import gammaln, polygamma
from scipy import special

import logging
import multiprocessing as mp
from itertools import repeat
import h5py
import tempfile

import os
import subprocess
from datetime import date

import re
import pdb
import time


"""rvd29.py: put Jeffrey's Prior on Mj."""



def main():
    import argparse

    # Populate our options, -h/--help is already there.
    description = """
                    RVD is a hierarchical bayesian model for identifying
                    rare variants from short-read sequence data. """
                    
    # create the top-level parser
    argp = argparse.ArgumentParser(prog='rvd', description=description)
    argp.add_argument('--version', action='version', version='%(prog)s 2.7')
    argp.add_argument('-v', '--verbose', dest='verbose', action='count',
                help="increase verbosity (specify multiple times for more)")
    
    # argp.add_argument('cmd', action='store', nargs='*',
    #             choices=['gen', 'gibbs'])
    subparsers = argp.add_subparsers(help='sub-command help')
    
    # create subparser for gibbs fitting
    argpGibbs = subparsers.add_parser('gibbs', 
                        help='fit the RVD model using Gibbs sampling')
    argpGibbs.add_argument('dcfile', nargs='+',
                        help='depth chart file name')
    argpGibbs.add_argument('-o', dest='outputFile', 
                default='output.hdf5',
                help='output HDF5 file name')
    argpGibbs.add_argument('-p', '--pool', type=int, default=None,
                help='number of workers in multithread pool')
    argpGibbs.set_defaults(func=gibbs)
                
                
    # create subparser to compare two model files
    argpTest = subparsers.add_parser('test', 
                        help='test if case error rate is greater than control by T')
    argpTest.add_argument('controlHDF5Name',
                help='control model file (HDF5)')
    argpTest.add_argument('caseHDF5Name',
                help='case model file (HDF5)')
    argpTest.add_argument('-T', type=float, default=0.005,
                help='threshold for computing variant probability (default=0.005)')
    argpTest.add_argument('-N', type=int, default=1000,
                help='Monte-Carlo sample size (default=1000)')
    argpTest.add_argument('-o', '--output', dest='outputFile', nargs='?', 
                default='test.hdf5')
    argpTest.set_defaults(func=test_main)
        
                
    # create subparser to sample the model
    argpGen = subparsers.add_parser('gen', 
                        help='sample data from the RVD model')
    argpGen.add_argument('input', nargs='+')
    argpGen.add_argument('-o', '--output', dest='outputFile', nargs='?', 
                default='output.hdf5')
                
    # Parse the arguments (defaults to parsing sys.argv)
    args = argp.parse_args()

    # TODO check what came in on the command line and call optp.error("Useful message") to exit if all is not well

    log_level = logging.WARNING  # default
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    
    # Set up basic configuration, out to stderr with a reasonable format
    logging.basicConfig(level=log_level,
                        format='%(levelname)s:%(module)s:%(message)s')
                        
    # Do actual work here
    args.func(args)
    
    
def gibbs(args):
    """ Top-level function to use gibbs sampling on a set of depth chart files
    """
    (r, n, loc, refb) = load_depth(args.dcfile)
    (phi, theta_s, mu_s, M_s) = mh_sample(r, n)
    save_model(args.outputfile, phi, mu=mu_s, M=M_s, theta=theta_s, r=r, n=n, loc=loc,
               refb=refb)

def test_main(args):
    test(args.controlHDF5Name, args.caseHDF5Name, args.T, args.N, args.outputFile)

def test(controlHDF5Name, caseHDF5Name, T=0.0005, N=1000, outputFile=None):
    """ Top-level function to test for variants.
    """
    
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    # Load the Case and Control Model files
    logging.debug(controlHDF5Name)
    (controlPhi, controlM, controlTheta, controlMu, controlLoc, controlR, controlN) = load_model(controlHDF5Name)
    (casePhi, caseM, caseTheta, caseMu, caseLoc, caseR, caseN) = load_model(caseHDF5Name)
    
    # Extract the common locations in case and control
    caseLocIdx = [i for i in xrange(len(caseLoc)) if caseLoc[i] in controlLoc]
    controlLocIdx = [i for i in xrange(len(controlLoc)) if controlLoc[i] in caseLoc]

    caseMu = caseMu[caseLocIdx,:]
    controlMu = controlMu[controlLocIdx,:]
    caseR = caseR[:,caseLocIdx,:]
    controlR = controlR[:,controlLocIdx,:]
    caseLoc = caseLoc[caseLocIdx]
    controlLoc = controlLoc[controlLocIdx]
    J = len(caseLoc)
    
    with h5py.File(controlHDF5Name, 'r') as f:
        refb = f['/refb'][...]
        f.close()
    refb = refb[controlLocIdx]
    
    
    # Sample from the posterior Z = muCase - muControl
    # Adjusting for baseline error rate for case and control
    (Z, caseMuS, controlMuS) = sample_post_diff(caseMu-casePhi['mu0'], controlMu-controlPhi['mu0'], N)
    
    # Posterior Prob that muCase is greater than muControl by T
    postP = bayes_test(Z, [(T, np.inf)]) 
    
    # chi2 test for goodness-of-fit to a uniform distribution for non-ref bases
    nRep = caseR.shape[0]
    chi2Prep = np.zeros((J,nRep))
    chi2P = np.zeros(J)
    for j in xrange(J):
        chi2Prep[j,:] = np.array([chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
        if np.any(np.isnan(chi2Prep[j,:])):
            chi2P[j] = np.nan
        else:
           chi2P[j] = 1-ss.chi2.cdf(-2*np.sum(np.log(chi2Prep[j,:] + np.finfo(float).eps)), 2*nRep) # combine p-values using Fisher's Method

    # Call variants that have postP > 0.95 and chi2pvalue < 0.05 after Bonferroni correction
    call = []
    altb = []
    for i in xrange(J):
        r = np.squeeze(caseR[:,i,:]) # replicates x bases
        
        # Make a list of the alternate bases for each replicate
        acgt_r = ['A','C','G','T']
        del acgt_r[ acgt[refb[i]] ]

        altb_r = [acgt_r[x] for x in np.argmax(r, axis=1)]
        
        if postP[i] >0.95 and chi2P[i] < 0.05/J: # Bonferroni Correction
##        if postP[i] >0.95: 
            altb.append(altb_r[0]) # TODO: find a better way to report all alternate bases
            call.append(True)
        else:
            altb.append(None)
            call.append(False)
     
    if outputFile is not None:
        # Save the test results
        with h5py.File(outputFile+'.hdf5', 'w') as f:
            f.create_dataset('loc', data=caseLoc)
            f.create_dataset('postP', data=postP)
            f.create_dataset('T', data=T)
            f.create_dataset('chi2pvalue',data=chi2P)
            f.close()
        write_vcf(outputFile+'.vcf', caseLoc, call, refb, altb, np.mean(caseMu, axis=1), np.mean(controlMu, axis=1))
    return caseLoc, caseMu, controlMu, postP, chi2P, call

def write_vcf(outputFile, loc, call, refb, altb, caseMu, controlMu):
    """ Write high confidence variant calls to VCF 4.2 file.
    """
    
    #TODO: get dbSNP id for chrom:pos
    J = len(loc)
    
    today=date.today()
    
    chrom = [x.split(':')[0][3:] for x in loc]
    pos = [int(x.split(':')[1]) for x in loc]
    vcfF = open(outputFile,'w')
    
    print("##fileformat=VCFv4.1", file=vcfF)
    print("##fileDate=%0.4d%0.2d%0.2d" % (today.year, today.month, today.day), file=vcfF)
    
    print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">", file=vcfF)
    print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">", file=vcfF)
    
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=vcfF)
    for i in xrange(J):
        if call[i]:
            print("%s\t%d\t.\t%c\t%s\t.\tPASS\tCOAF=%0.3f;CAAF=%0.3f" % (chrom[i], pos[i], refb[i], altb[i], controlMu[i]*100.0, caseMu[i]*100.0), file=vcfF)
   
    vcfF.close()
    
'''
def sample_run():
    n = 1000
    J = 10
    phi = {'mu0': 0.20, 'M0': 2e3, 'a': 1e6, 'b': 1}
    r, theta, mu, M = generate_sample(phi, n=n, J=J, seedint=10)
    r[:, int(J / 2)] = n * np.array([0.50, 0.55, 0.45])

    phi, theta_s, mu_s = mh_sample(r, n, 
                                        nsample=100, 
                                        thin=0, 
                                        burnin=0)							
'''
						
def sample_run():
    n = 100
    J = 100
    M = [100 for i in range(J)]
    r, theta, mu, M = generate_sample(M, J, n);
    r[:, int(J/2)] = n * np.array([0.50, 0.55, 0.45])

    phi, theta_s, mu_s = mh_sample(r, n, 
                                        nsample=100, 
                                        thin=0, 
                                        burnin=0)


def load_model(h5Filename):
    """ Returns the RVD29 model samples and parameters.
    Takes an hdf5 filename and returns phi and other parameters
    """

    out = []

    with h5py.File(h5Filename, 'r') as h5file:
        # Load phi - it always exists
        phi = {'mu0': h5file['phi/mu0'][()],
               'M0': h5file['phi/M0'][...]}
        out.append(phi)
        
        # Load M if it exists
        if u"M" in h5file.keys():
            M = h5file['M'][...]
            out.append(M)
            
        # Load theta if it exists
        if u"theta" in h5file.keys():
            theta = h5file['theta'][...]
            out.append(theta)
            
        # Load mu if it exists
        if u"mu" in h5file.keys():
            mu = h5file['mu'][...]
            out.append(mu)
            
        # Load loc if it exists
        if u"loc" in h5file.keys():
            loc = h5file['loc'][...]
            out.append(loc)
            
        # Load r if it exists
        if u"r" in h5file.keys():
            r = h5file['r'][...]
            out.append(r)

        if u"n" in h5file.keys():
            n = h5file['n'][...]
            out.append(n)

            
    return tuple(out)
      
def save_model(h5Filename, phi, mu=None, M=None, theta=None, r=None, n=None, loc=None, refb=None):
    """ Save the RVD29 model samples and parameters """
    
    # TODO add attributes to hdf5 file
    h5file = h5py.File(h5Filename, 'w')
    
    # Save the model parameters (phi)
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', data=phi['mu0'])
    h5file['phi'].create_dataset('M0', data=phi['M0'])
    h5file.create_dataset('M', data=M, 
                                      chunks=True, 
                                      fletcher32=True, 
                                      compression='gzip')
    
    # Save the latent variables if available.
    if mu is not None:
        h5file.create_dataset('mu', data=mu, 
                              chunks=True, fletcher32=True, compression='gzip')
    if theta is not None:
        h5file.create_dataset('theta', data=theta, 
                              chunks=True, fletcher32=True, compression='gzip')
    
    # Save the data used for fitting the model if available
    if r is not None:
        h5file.create_dataset('r', data=r, 
                              chunks=True, fletcher32=True, compression='gzip')
    if n is not None:
        h5file.create_dataset('n', data=n, 
                              chunks=True, fletcher32=True, compression='gzip')

    # Save the reference data
    if loc is not None:
        h5file.create_dataset('loc', data=loc, 
                              chunks=True, fletcher32=True, compression='gzip')
    if refb is not None:
        h5file.create_dataset('refb', data=refb)

    h5file.close()


def generate_sample(mu0=0.25, M0=10, n=100, N=3, J=100, seedint=None):
    """Returns a sample with reads, N replicates, and
    J locations. The parameters of the model are M0, mu0. In addition,
	M is initialed randomly."""
    if seedint is not None: 
        np.random.seed(seedint)
        
    # Draw J location-specific error rates from a Beta
    #mu0 = float(mu0)
    alpha0 = M0*mu0
    beta0 = M0*(1-mu0)
    mu = ss.beta.rvs(alpha0, beta0, size=J)
    
    #generate initial sample M
    M = np.dot(1000,np.random.rand(J))
    #print (M)
    
    # Draw J location-specific precisions from Jeffreys' Prior
    for j in xrange(0, J):
        M[j] = polygamma(1, mu[j]*M[j])*(mu[j]**2) + polygamma(1, (1-mu[j]*M[j]))*(1-mu[j])**2 -polygamma(1,M[j])
        M[j] = np.sqrt(M[j])
    #print (M)
    
    # Draw sample error rate and error count
    theta=np.zeros((N,J))
    r = np.zeros((N,J))
    for j in xrange(0, J):
        alpha = mu[j]*M[j]
        beta = (1-mu[j])*M[j]
        theta[:,j] = ss.beta.rvs(alpha, beta, size=N)
        r[:,j] = ss.binom.rvs(n, theta[:,j])
    return r, theta, mu, M


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

    theta = r/(n + np.finfo(np.float).eps) # make sure this is non-truncating division
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps)
    
    M = (mu*(1-mu))/(np.var(mu, 0) + np.finfo(np.float).eps )
    
    phi = {'mu0':mu0, 'M0':M0}
    return phi, mu, theta, M
    
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
    
def sampleMuMH(theta, mu0, M0, M, mu=ss.beta.rvs(1, 1), burnin=0, mh_nsample=1, thin=0, pool=None):
    """ Return a sample of mu with parameters mu0 and M0.
    """
    if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
    elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
        
    alpha0 = mu0*M0 + np.finfo(np.float).eps
    beta0 = (1-mu0)*M0 + np.finfo(np.float).eps
    
    Qsd = mu/10
    mu_s = np.zeros( (mh_nsample, J) ) 
    for ns in xrange(0, mh_nsample):
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
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*mh_nsample):], 0)
    if thin > 0:
        mu_s = np.delete(mu_s, np.s_[::thin], 0)
    
    return mu_s
    

def sampleLocMMH(args):
    mu, Qsd, theta, M = args
    # Sample from the proposal distribution
    while True:
        M_p = ss.norm.rvs(M, Qsd)
        if M_p > 0: break
        
    # Log-likelihood for the proposal M
    alpha_p = mu*M_p + np.finfo(np.float).eps
    beta_p = (1-mu)*M_p + np.finfo(np.float).eps
    logPM_p = np.sum(beta_log_pdf(theta, alpha_p, beta_p)) \
                + np.log(( np.sqrt(polygamma(1, np.dot(mu,M_p))*(mu**2) + polygamma(1, np.dot(1-mu,M_p))*(1-mu)**2 -polygamma(1,M_p))))
                    
    # Log-likelihood for the current M
    alpha = mu*M + np.finfo(np.float).eps
    beta = (1-mu)*M + np.finfo(np.float).eps
    logPM = np.sum(beta_log_pdf(theta, alpha, beta)) \
                + np.log(( np.sqrt(polygamma(1, np.dot(mu,M))*(mu**2) + polygamma(1, np.dot(1-mu,M))*(1-mu)**2 -polygamma(1,M))))

    # Accept new M if it increases posterior pdf or by probability
    loga = logPM_p - logPM
    # if j==0: print (M_p, M[j], logPM_p, logPM)
    if (loga > 0 or np.log(np.random.random()) < loga): 
        M = np.copy(M_p)
    return M

def sampleMMH(M, theta, mu, burnin=0, mh_nsample=1, thin=0, pool=None):
    """ Return a sample of M with parameters theta, mu, phi['M']"""
    if np.ndim(theta) == 1: (N, J) = (1, np.shape(theta)[0])
    elif np.ndim(theta) > 1: (N, J) = np.shape(theta)
        
    
    Qsd = M/10
    M_s = np.zeros( (mh_nsample, J) ) 
    for ns in xrange(0, mh_nsample):
        if pool is not None:
            args = zip(mu, Qsd, theta.T, M)
            M = pool.map(sampleLocMMH, args)
        else:
            for j in xrange(0, J):
                args = (mu[j], Qsd[j], theta[:,j], M[j])
                M[j] = sampleLocMMH(args)
        # Save the new sample
        M_s[ns, :] = np.copy(M)

    if burnin > 0.0:
        M_s = np.delete(M_s, np.s_[0:np.int(burnin*mh_nsample):], 0)
    if thin > 0:
        M_s = np.delete(M_s, np.s_[::thin], 0)
    
    return M_s



def mh_sample(r, n, gibbs_nsample=10000, mh_nsample=10, burnin=0.2, thin=2, pool=None):
    
    """Return MAP parameter and latent variable estimates obtained by 
    Metropolis-Hastings sampling.
    By default, sample 10000 M-H with a 20% burn-in and thinning factor of 2. 
    Stop when the change in complete data log-likelihood is less than 0.01%.
    """
    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    elif np.ndim(r) == 3: 
        r = np.sum(r, 2) # sum over non-reference bases
        N, J = np.shape(r)
    
    # Initialize a hdf5 file for logging model progress
    h5Filename = tempfile.NamedTemporaryFile(suffix='.hdf5').name
    logging.debug("Storing temp data in %s" % h5Filename)
    h5file = h5py.File(h5Filename, 'w')
    h5file.create_group('phi')
    h5file['phi'].create_dataset('mu0', (1,), dtype='f')
    h5file['phi'].create_dataset('M0', (1,), dtype='f')
    h5file.create_dataset('theta_s', (N, J, gibbs_nsample), dtype='f')
    h5file.create_dataset('mu_s', (J, gibbs_nsample), dtype='f')
    h5file.create_dataset('M_s', (J, gibbs_nsample), dtype='f')

    # Initialize estimates using MoM
    phi, mu, theta, M = estimate_mom(r, n)
    logging.debug("MoM: mu0 = %0.3e; M0 = %0.3e." % (phi['mu0'], phi['M0']) )
    #logging.debug("MoM: a = %0.33; b = %0.3e." % (phi['a'], phi['b']))

    # Correct MoM estimates to be non-trivial
    mu[mu < np.finfo(np.float).eps*1e4] = phi['mu0']
    theta[theta < np.finfo(np.float).eps*1e4] = phi['mu0']
    M[M < np.finfo(np.float).eps *1e4] = 1
    
    # Sample theta, mu by gibbs sampling
    theta_s = np.zeros( (N, J, gibbs_nsample) )
    mu_s = np.zeros( (J, gibbs_nsample) )
    M_s = np.zeros( (J, gibbs_nsample) )
    
    for i in xrange(gibbs_nsample):
        if i % 100 == 0 and i > 0: logging.debug("Gibbs Iteration %d" % i)
            
        # Draw samples from p(theta | r, mu, M) by Gibbs
        alpha = r + mu*M +  + np.finfo(np.float).eps
        beta = (n - r) + (1-mu)*M + np.finfo(np.float).eps
        theta = ss.beta.rvs(alpha, beta)
        
        # Draw samples from p(mu | theta, mu0, M0) by Metropolis-Hastings
        mu_mh = sampleMuMH(theta, phi['mu0'], phi['M0'], M, mu=mu, mh_nsample=mh_nsample, pool=pool)
        
        
        # Draw sample from p(M | theta, a, b) by Metropolis-Hastings
        M_mh = sampleMMH(M, theta, mu=mu, mh_nsample=mh_nsample, pool=pool)

        mu = np.median(mu_mh, axis=0) 
        M = np.median(M_mh, axis=0)
        
        
        # Store the sample
        theta_s[:,:,i] = np.copy(theta)
        mu_s[:,i] = np.copy(mu)
        M_s[:,i] = np.copy(M)
        # Update parameter estimates
        # phi['mu0'] = np.mean(mu)
        # phi['M0'] = (phi['mu0']*(1-phi['mu0']))/(np.var(mu) + np.finfo(np.float).eps)
        # TODO update for M
        
        # Store the current model
        h5file['phi']['mu0'][0] = phi['mu0']
        h5file['phi']['M0'][0] = phi['M0']
        h5file['M_s'][:,i] = M
        h5file['theta_s'][:,:,i] = theta
        h5file['mu_s'][:,i] = mu
        h5file.flush()
    
    # Apply the burn-in and thinning
    if burnin > 0.0:
        mu_s = np.delete(mu_s, np.s_[0:np.int(burnin*gibbs_nsample):], 1)
        theta_s = np.delete(theta_s, np.s_[0:np.int(burnin*gibbs_nsample):], 2)
        # identity the paramter 
        M_s = np.delete(M_s, np.s_[0:np.int(burnin*gibbs_nsample):], 1)
    if thin > 1:
        mu_s = np.delete(mu_s, np.s_[::thin], 1)
        theta_s = np.delete(theta_s, np.s_[::thin], 2)
        # identity the paramter 
        M_s = np.delete(M_s, np.s_[::thin], 1)
    
    h5file.close()
    return (phi, theta_s, mu_s, M_s)


def beta_log_pdf(x, a, b):
    return gammaln(a+b) - gammaln(a) - gammaln(b) \
            + (a-1)*np.log(x+np.finfo(np.float).eps) \
            + (b-1)*np.log(1-x)


def ll(phi, r):
    """ Return the log-likelihood of the data, r, under the model, phi.
    """
    pass

def make_pileup(bamFileName, fastaFileName, region):
    """ Creates a pileup file using samtools mpileup in /pileup directory.
    """
    
    # Check that BAM file exists
    assert os.path.isfile(bamFileName), "BAM file does not exist: %s" % bamFileName
    
    # Check that FASTA reference file exists
    assert os.path.isfile(fastaFileName), "FASTA file does not exist: %s" % fastaFileName
    
    # Create pileup directory if it doesn't exist
    if not os.path.isdir("pileup"):
        os.makedirs("pileup")
    
    # Format the samtools call
    callString = ["samtools", "mpileup", "-d", "1000000", "-r", "%s" % region,
                  "-f", "%s" % fastaFileName, "%s" % bamFileName]
                  
    # Remove the extension from the bam filename and replace with .pileup
    pileupFileName = bamFileName.split("/")[-1]
    pileupFileName = os.path.join("pileup", 
                                  "%s.pileup" % pileupFileName.split(".", 1)[0])
    
    # Run samtools pileup only if the file doesn't already exist.
    #try:
    #    with open(pileupFileName, 'r'):
    #        logging.debug("Pileup file exists: %s" % pileupFileName)
    #except IOError:
    logging.debug("[call] %s", " ".join(callString))
    with open(pileupFileName, 'w') as fout:
        subprocess.call(callString, stdout=fout)
    return pileupFileName

def make_depth(pileupFileName):
    """ Generates a depth chart file for each pileup file and stores it in the
        /depth_chart directory. The folder will be created if it doesn't exist.
    """
    
    if not os.path.isdir("depth_chart"):
        os.makedirs("depth_chart")
    
    # TODO replace this with a python version.
    callString = ["../../bin/pileup2dc", "%s" % pileupFileName]

    dcFileName = pileupFileName.split("/")[-1]
    dcFileName = os.path.join("depth_chart", 
                                  "%s.dc" % dcFileName.split(".", 1)[0])
    #try:
    #    with open(dcFileName, 'r'):
#		logging.debug("Depth chart file exists: %s" % dcFileName) 
    #except IOError:
    logging.debug("Converting %s to depth chart." % pileupFileName)
    with open(dcFileName, 'w') as fout:
        subprocess.call(callString, stdout=fout)
    return dcFileName

def load_depth(dcFileNameList):
    """ Return (r, n, location, reference base) for a list of depth charts. The
        variable r is the error read depth and n is the total read depth.
    """
    r=[]; n=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    loc = []
    refb = {}
    cd = []
    for dcFileName in dcFileNameList:
        with open(dcFileName, 'r') as dcFile:
            header = dcFile.readline().strip()
            dc = dcFile.readlines()
            dc = [x.strip().split("\t") for x in dc]
            
            loc1 = [x[1]+':'+str(x[2]).strip('\000') for x in dc if x[4] in acgt.keys()]
	    loc.append( loc1 )
            
            refb1 = dict(zip(loc1, [x[4] for x in dc if x[4] in acgt.keys()]))
            refb.update(refb1)
            cd.append( dict(zip(loc1, [map(int, x[5:9]) for x in dc if x[4] in acgt.keys()])) )
            
    loc = list(reduce(set.intersection, map(set, loc)))

    def stringSplitByNumbers(x):
        r = re.compile('(\d+)')
        l = r.split(x)
        return [int(y) if y.isdigit() else y for y in l]

    loc = sorted(loc,key = stringSplitByNumbers)
    logging.debug(loc)
    refb = [refb[k] for k in loc]
    
    J = len(loc)
    N = len(dcFileNameList)
    for i in xrange(0, N):
        logging.debug("Processing %s" % dcFileNameList[i])
        c = np.array( [cd[i][k] for k in loc] )
        n1 = np.sum(c, 1)
        #r1 = np.zeros(J)
        refIdx=np.zeros(J)

        for j in xrange(0,J):
            #r1[j] = n1[j] - c[j, acgt[refb[j]]]
            refIdx[j] = 4*j+acgt[refb[j]]
        c = np.delete(c, refIdx, None)
        c = np.reshape(c, (J, 3) )
        #r.append(r1)
        n.append(n1)
        r.append(c)
    r = np.array(r)
    n = np.array(n)

    return (r, n, loc, refb)

def chi2test(X, lamda=2.0/3, pvector=np.array([1.0/3]*3)):
    """ Do chi2 test to decide how well the error reads fits uniform multinomial distribution. P-value returned.
        lamda=1 Pearson's chi-square
        lamda=0 the log likelihood ratio statistic/ G^2
        lamda=-1/2 Freeman-Tukey's F^2
        lamda=-1  Neyman modified chi-square
        lamda=-2  modified G^2
    """
    X=np.array(X)

    nsum=np.sum(X)
    if nsum == 0: return np.nan # return NaN if there are no counts
    E=nsum*pvector


    if lamda==0 or lamda==-1:
        C=2.0*np.sum(X*np.log(X*1.0/E))
    else:
        C=2.0/(lamda*(lamda+1))*np.sum(X*((X*1.0/E)**lamda-1))
        
    df=len(pvector)-1
    #p=scipy.special.gammainc(C,df)
    # p=1-gammainc(df/2,C/2)
    p = 1 - ss.chi2.cdf(C, df) 
    return(p)
 
def sample_post_diff(muCaseG, muControlG, N):
    """ Return N samples from the posterior distribution for 
         u_j|r_case - u_j|r_control. """
    
    nCase = muCaseG.shape[1]
    nControl = muControlG.shape[1]
    
    caseSample = np.random.choice(nCase, size=N, replace=True)
    controlSample = np.random.choice(nControl, size=N, replace=True)
    
    muCaseS = muCaseG[:, caseSample]
    muControlS = muControlG[:, controlSample]

    Z = muCaseS - muControlS

    return (Z, muCaseS, muControlS)
 
def bayes_test(Z, roi):
    """ Return posterior probabilities in regions defined in list of tuples (roi)
        from samples in columns of Z. """

    (J,N)=np.shape(Z)
    
    nTest = len(roi) # get the number of regions to compute probabilities 
    
    p = np.zeros((J,nTest))
    for i in xrange(nTest):
        for j in xrange(J):
		p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] >= roi[i][0]), (Z[j,:] < roi[i][1]) ) ) ) / N

    return p
    
if __name__ == '__main__':
    main()
    
