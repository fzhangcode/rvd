# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import numpy as np

import scipy.stats as ss
import scipy.optimize as so
from scipy.special import gammaln, psi, betaln
from scipy import linalg, integrate

#import pandas as pd

import multiprocessing as mp
from itertools import repeat

import h5py as h5
import tempfile
import logging
import time
from datetime import datetime

import pdb

def main():
    log_level = logging.DEBUG # default logging level
    logging.basicConfig(level=log_level, format='%(levelname)s:%(module)s:%(message)s')

    ## Generate simulation data

    ## model optimization  
    # (phi, q) = varlap_opt(y, K=2, seed=10241978, pool=pool)

    ## save the parameters.
    # save_model('model.hdf5', y, phi, q)
    

    

## compute sufficient statistics
def EqlogTheta(delta):
	return psi(delta[0]) - psi(np.sum(delta))

def Eqlog1_Theta(delta):
	return psi(delta[1]) - psi(np.sum(delta))

def EqMu(gam):
	return gam[0] / (np.sum(gam)) # eps?

def EqlogMu(gam):
    return psi(gam[0]) - psi(np.sum(gam))

def Eqlog1_Mu(gam):
	return psi(gam[1]) - psi(np.sum(gam))

# def EqlogGamma(gam, M):
# 	# Eqlog\gam(muM)
# 	N = 7
# 	mu = range(0, 1 , 2^N+1)
# 	y = [kernel(x, gam , M) for x in mu]

# 	return integrate.romb(y, dx = 0.1/(2^N+1))

def EqlogGamma(gam, M):
    logGamma = integrate.quad(kernel, 0, 1, args = (gam, M))
    return logGamma[0]

def kernel(mu, gam, M):
    return ss.beta.pdf(mu, gam[0], gam[1])*betaln(mu*M, (1-mu)*M)

## compute entropy
def BetaEntropy(x):
	# To compute EqlogQmu and EqlogQtheta
	return betaln(x[0],x[1]) - (x[0]-1) * psi(x[0]) - (x[1] - 1) * psi(x[1]) +(x[0] + x[1] -2) * psi(x[0] + x[1])


## compute ELBO
def ELBO(r, n, M, mu0, M0, delta, gam):
    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    
    Mu = np.array([EqMu(gam[j,:]) for j in xrange(J)])
    logMu = np.array([EqlogMu(gam[j,:]) for j in xrange(J)])
    log1_Mu = np.array([Eqlog1_Mu(gam[j,:]) for j in xrange(J)])

    logTheta = np.zeros((N,J))
    log1_Theta = np.zeros((N,J))

    for j in xrange(J):
        for i in xrange(N):
            logTheta[i,j] = EqlogTheta(delta[i,j,:])
            log1_Theta[i,j] = Eqlog1_Theta(delta[i,j,:])


    EqlogPr = 0.0
    for j in xrange(J):
        for i in xrange(N):
            EqlogPr += - betaln(r[i,j] + 1, n[i,j] - r[i,j] +1)
            EqlogPr += r[i,j]*logTheta[i,j] + (n[i,j] - r[i,j]) * logTheta[i,j]

    EqlogPtheta = 0.0
    for j in xrange(J):
        
        EqlogPtheta += N*EqlogGamma(gam[j,:], M[j])
        for i in xrange(N):
            EqlogPtheta += (M[j]* Mu[j]- 1)*logTheta[i,j] +\
            (M[j]*(1 - Mu[j]) - 1)*log1_Theta[i,j]

    EqlogPmu = -J * betaln(mu0*M0, (1-mu0)*M0)
    for j in xrange(J):
        EqlogPmu += (M0*mu0-1)*logMu[j] + (M0*(1-mu0)-1)*log1_Mu[j]

    EqlogQtheta = 0.0
    for j in xrange(J):
        for i in xrange(N):
            EqlogQtheta += BetaEntropy(delta[i,j,:])

    EqlogQmu = 0.0
    for j in xrange(J):
        EqlogQmu += BetaEntropy(gam[j,:])

    return EqlogPr + EqlogPtheta + EqlogPmu - EqlogQtheta - EqlogQmu

def ELBO_delta_ij(r, n, M, delta, gam):
    ## partial ELBO from replicate i position j
    ## ELBO used to optimize delta
    ## Commented out all items that don't depend on delta
  
    Mu = EqMu(gam)
    logTheta = EqlogTheta(delta)
    log1_Theta = Eqlog1_Theta(delta)

    EqlogPr = betaln(r+1, n-r+1) + r*logTheta + (n - r)*logTheta

    EqlogPtheta = EqlogGamma(gam, M) + (M*Mu - 1)*logTheta + (M*(1-Mu)-1)*log1_Theta

    EqlogQtheta = BetaEntropy(delta)

    return EqlogPr + EqlogPtheta  - EqlogQtheta

def neg_ELBO_delta_ij(logdelta, gam, r, n, M):
    return -ELBO_delta_ij(r, n, M, np.exp(logdelta), gam)

def opt_delta_ij(args):
    r, n, M, delta, gam = args

    bnds = [[-7, 7]]*2# limit delta to [0.0001, 10000]
    args=(gam, r, n, M)

    logdelta = opt_par(neg_ELBO_delta_ij, np.log(delta), args, bnds, 'delta')
    delta = np.exp(logdelta)
    return delta

def opt_delta(r, n, M, delta, gam, pool = None):
    logging.debug("Optimizing delta")

    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)

    st = time.time()
    if pool is not None:
        for i in xrange(N):
            args = zip (r[i,:], n[i,:], M, delta[i,:,:], gam)
            temp = pool.map(opt_delta_ij, args)
            delta = np.array(temp)

    else:
        logging.debug('Optimizing delta in single thread')
        for i in xrange(N):
            for j in xrange(J):
                args = (r[i,j], n[i,j], M[j], delta[i,j,:], gam[j,:])
                delta[i,j,:] = opt_delta_ij(args)

def ELBO_gam_j( M, mu0, M0, delta, gam):
    ## partial ELBO depending on gam from each position j
    ## ELBO used to gam

    if np.ndim(delta) == 1: N = 1
    elif np.ndim(delta) == 2: N= np.shape(delta)[0]

    Mu = EqMu(gam)
    logMu = EqlogMu(gam)
    log1_Mu = Eqlog1_Mu(gam)

    logTheta = np.zeros((N,1))
    log1_Theta = np.zeros((N,1))

    for i in xrange(N):
        logTheta[i] = EqlogTheta(delta[i,:])
        log1_Theta[i] = Eqlog1_Theta(delta[i,:])

    EqlogPtheta = N*EqlogGamma(gam,M)
    for i in xrange(N):
        EqlogPtheta += (M*M-1) * logTheta[i] + (M*(1-Mu)-1)*log1_Theta[i]

    EqlogPmu= -betaln(mu0*M0, (1-mu0)*M0)+ (M0*mu0-1)*logMu + (M0*(1-mu0)-1)*log1_Mu

    EqlogQmu = BetaEntropy(gam)

    return  EqlogPtheta + EqlogPmu - EqlogQmu

def neg_ELBO_gam_j(loggam, delta, M, mu0, M0):
    return -ELBO_gam_j(M, mu0, M0, delta, np.exp(loggam))

def opt_gam_j(args):
    M, mu0, M0, delta, gam = args
    bnds = [[-7, 7]]*2# limit delta to [0.0001, 10000]
    args = (delta, M, mu0, M0)
    loggam = opt_par(neg_ELBO_gam_j, np.log(gam), args, bnds, 'gamma')
    gam = np.exp(loggam)
    return gam

def opt_gam(M, mu0, M0, delta, gam, pool = None):
    logging.debug("Optimizing gam")

    if np.ndim(gam) == 1: J=1
    elif np.ndim(gam) == 2: J=np.shape(gam)[0]

    st = time.time()

    if pool is not None:
        args = zip( M, repeat(mu0,J), repeat(M0,J), np.transpose(delta,axes=(1,0,2)),gam)
        gam = pool.map(opt_gam_j, args)

    else:
        for j in xrange(J):
            args = ( M[j], mu0, M0, delta[:,j,:], gam[j] )
            gam[j] = opt_gam_j(args)

    logging.debug('Gamma update elapsed time is %0.3f sec for %d samples.' % (time.time() - st, J))
    return gam

def ELBO_0(mu0, M0, gam):
    ## Items in ELBO depends on mu0 and M0
    ## FOr optimization of mu0 and M0

    J = gam.shape[0]
    
    logMu = np.array([EqlogMu(gam[j,:]) for j in xrange(J)])

    log1_Mu = np.array([Eqlog1_Mu(gam[j,:]) for j in xrange(J)])

    EqlogPmu = -J * betaln(mu0*M0, (1-mu0)*M0)
    for j in xrange(J):
        EqlogPmu += (M0*mu0-1)*logMu[j] + (M0*(1-mu0)-1)*log1_Mu[j]

    return EqlogPmu

def neg_ELBO_mu0(mu0, M0, gam):
    return -ELBO_0(mu0, M0, gam)

def opt_mu0(mu0, M0, gam):
    logging.debug("Optimizing mu0")
    bnds = np.array([[0.01,0.99]])
    args=(M0, gam)
    mu0 = opt_par(neg_ELBO_mu0, mu0, args, bnds, 'mu0' )

    return mu0

def neg_ELBO_M0(logM0, mu0, gam):
    return -ELBO_0(mu0, np.exp(logM0), gam)
    
def opt_M0(mu0, M0, gam):
    logging.debug("Optimizing M0")

    bnds = np.array([[-7,7]])
    args = (mu0, gam)
    logM0 = opt_par(neg_ELBO_M0, np.log(M0), args, bnds, 'M0' )
    M0 = np.exp(logM0)
    return M0

def ELBO_M_j(M, delta, gam):
    ## partial ELBO depending on M from each position j 
    ## ELBO used to optimize M 

    if np.ndim(delta) == 1: N = 1
    elif np.ndim(delta) == 2: N= np.shape(delta)[0]

    Mu = EqMu(gam)

    logTheta = np.zeros((N,1))
    log1_Theta = np.zeros((N,1))

    for i in xrange(N):
        logTheta[i] = EqlogTheta(delta[i,:])
        log1_Theta[i] = Eqlog1_Theta(delta[i,:])

    EqlogPtheta = N*EqlogGamma(gam,M)
    for i in xrange(N):
        EqlogPtheta += (M*M-1) * logTheta[i] + (M*(1-Mu)-1)*log1_Theta[i]

    return  EqlogPtheta 

def neg_ELBO_M_j(logM, delta, gam):
    return -ELBO_M_j(np.exp(logM), delta, gam)

def opt_M_j(args):

    (M, delta, gam) = args
    bnds = np.array([[-1, 11]]) # limit delta to [0.0001, 10000]
    args = (delta, gam)

    logM = opt_par(neg_ELBO_M_j, np.log(M), args, bnds, 'M')
    M = np.exp(logM)

    return M

def opt_M(M, delta, gam, pool = None):
    logging.debug("Optimizing M")

    if np.ndim(M) == 1: J = 1
    elif np.ndim(M) == 2: J= np.shape(M)[0]

    if pool is not None:
        args = zip(M, np.transpose(delta, axes =(1,0,2)), gam)
        M = pool.map(opt_M_j, args)

    else:
        for j in xrange(J):
            args = (M[j],delta[:,j,:],gam[j,:])

    return M

def opt_par(func, x, args, bnds, parlable):
    res = so.minimize(func, x, 
        args=args, bounds=bnds, 
        method='L-BFGS-B')

    # logging.debug("Inside of optimize function. got res")
    if res.success == False:
        res = so.minimize(func, x, 
            args=args, method='Nelder-Mead')

    if res.success == False:
        so.minimize(func, x, \
            args=args, method='BFGS')
    
    if res.success == False or np.any ( np.isnan(res.x) ):
        logging.warning("Could not optimize %s or %s is NaN." %(parlable, parlable))
        if np.ndim(bnds) == 1:
            size = 1
        else:
            size = np.shape(bnds)  
        x = np.random.uniform(low=np.amin(bnds), high=np.amax(bnds), size = size)      
        return x

    x = res.x
    return x

def ELBO_opt(r, n, phi = None, q = None, seed = None, pool = None):
    t = str(datetime.now)
    f = open('%s.txt' %t,'w')

    print("ELBO optimization trace starting from %s: \n" %t, file=f)

    if np.ndim(r) == 1: N, J = (1, np.shape(r)[0])
    elif np.ndim(r) == 2: N, J = np.shape(r)
    elif np.ndim(r) == 3: 
        r = np.sum(r, 2) 
        (N, J) = r.shape# sum over non-reference bases
    # r = r.T
    # n = n.T

    if seed is not None: np.random.seed(seed = seed)

    h5file = tempfile.NamedTemporaryFile(suffix='.hdf5')
    logging.info('Storing model updates in %s' % h5file.name)   

    ## Define optimization stopping criterion
    # MAXITER = 20
    # ELBOTOLPCT = 0.1    
    # MAXVARITER = 10
    # NORMTOL = 0.1
    MAXITER = 2
    ELBOTOLPCT = 0.1    
    MAXVARITER = 2
    NORMTOL = 0.1


    ## Initialize model parameters
    if phi is None:
        phi, mu, theta = estimate_mom(r, n)
    else:
        _, mu, theta = estimate_mom(r, n)
    mu0 = phi['mu0']
    M0 = phi['M0']
    M = phi['M']

    ## Initialize the variational parameters
   
    if q is None:
        # delta = np.transpose(np.array([r + mu*M+ np.finfo(np.float).eps, 
        #     (n - r) + (1-mu)*M+ np.finfo(np.float).eps]), axes = (1,2,0))
        delta = np.random.uniform(low = 0.1, high = 100, size = (N,J,2))
        gam = np.random.uniform(low=0.1, high=100, size = (J,2))
    else:
        delta = q['delta']
        gam = q['gam']

    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    q = {'delta':delta, 'gam':gam}
    #save_model('initial_value.hdf5', r, n, phi, q)  

    ## Initialize ELBO
    elbo = [ELBO(r, n, M, mu0, M0, delta, gam)]
    logging.info("Initial ELBO: %0.2f" % elbo[-1])


    ## Optimization
    moditer = 0
    delta_elbo_pct = np.inf

    while moditer < MAXITER and delta_elbo_pct > ELBOTOLPCT:
        # E-step: Update the variational distribution
        variter = 0
        var_elbo = [ elbo[-1] ]
        (norm_delta_delta, norm_delta_gam) = (np.inf, np.inf)

        while variter < MAXVARITER \
            and delta_elbo_pct > ELBOTOLPCT \
            and (norm_delta_delta > NORMTOL or norm_delta_gam > NORMTOL):

            #Store the previous parameter values
            (delta_prev, gam_prev) = (np.copy(delta), np.copy(gam))

            #Update the variational distribution
            gam = opt_gam( M, mu0, M0, delta, gam, pool = pool)
            delta = opt_delta(r, n, M, delta, gam, pool = pool)
            
            #Test for convergence
            var_elbo.append(ELBO(r, n, M, mu0, M0, delta, gam))
            delta_varelbo_pct = 100(var_elbo[-1] - var_elbo[-2])/abs(var_elbo[2])
            logging.info("Variational Step ELBO: %0.2f; Percent Change: %0.3f%%" % (var_elbo[-1], delta_varelbo_pct))
          
            norm_delta_delta = linalg.norm(delta - delta_prev)
            norm_delta_gam = linalg.norm(gam - gam_prev)
            logging.debug("||delta - delta_prev|| = %0.2f; ||gam - gam_prev|| = %0.2f" 
                % (norm_delta_delta, delta_norm_gam))

            variter += 1

        # M-step: Update model parameters
        mu0 = opt_mu0(mu0, M0, gam)
        M0 = opt_M0(mu0, M0, gam)
        M = opt_M(M, delta, gam, pool = pool)

        elbo.append(ELBO(r, n, M, mu0, M0, delta, gam))
        delta_elbo_pct = 100*(elbo[-1] - elbo[-2])/abs(elbo[-2])
        moditer += 1

        # ibic

        # Display results for debugging
        logging.info("Iteration %d of %d." % (moditer, MAXITER))
        logging.info("ELBO: %0.2f; Percent Change: %0.2f%%" \
                    % (elbo[-1], delta_elbo_pct))
        
        print("ELBO: %0.2f; Percent Change: %0.2f%%" % (elbo[-1], delta_elbo_pct), file = f)

        logging.info("M0 = %0.2e" % M0)
        logging.info("mu0 = %0.2f" % mu0)

        # Store the model for viewing
        phi = {'mu0':mu0, 'M0':M0, 'M':M}
        q = {'delta':delta, 'gam':gam}
        save_model(h5file.name, r, n, phi, q)

        phi = {'mu0':mu0, 'M0':M0, 'M':M}
        q = {'delta':delta, 'gam':gam}

    f.close()
    return(phi, q)

def estimate_mom(r, n):
    """ Return model parameter estimates using method-of-moments.
    """
    theta = r/(n + np.finfo(np.float).eps) # make sure this is non-truncating division
    if np.ndim(r) == 1: mu = theta
    elif np.ndim(r) > 1: mu = np.mean(theta, 0)
    
    mu0 = np.mean(mu)
    M0 = (mu0*(1-mu0))/(np.var(mu) + np.finfo(np.float).eps) + np.finfo(np.float).eps

    # estimate M. If there is only one replicate, set M as 10 times of M0.
    # If there is multiple replicates, set M according to the moments of beta distribution

    if np.shape(theta)[0] is 1:
        M = 10*M0*np.ones_like(mu)
    else:
        M = (mu*(1-mu))/(np.var(theta, 0) + np.finfo(np.float).eps ) 

    J = len(M)
    for i in xrange(J):
        if M[i] < 1:
            M[i] = 1
   


    phi = {'mu0':mu0, 'M0':M0, 'M':M}
    return phi, mu, theta

def save_model(fname, r, n, phi, q):

    f = h5.File(fname, 'w')
    
    f.create_dataset('r', data=r)
    f.create_dataset('n', data=n)

    f.create_group('phi')
    f['phi'].create_dataset('mu0', data=phi['mu0'])
    f['phi'].create_dataset('M0', data=phi['M0'])
    f['phi'].create_dataset('M', data=phi['M'])
    
    f.create_group('q')
    f['q'].create_dataset('delta', data=q['delta'])
    f['q'].create_dataset('gam', data=q['gam'])

    f.close()

def load_model(fname):
    
    f = h5.File(fname, 'r')
    
    r = f['r'][...]
    n = f['n'][...]

    phi = {}    
    phi['mu0'] = f['phi/mu0'][...]
    phi['M0'] = f['phi/M0'][...]
    phi['M'] = f['phi/M'][...]

    q = {}
    q['delta'] = f['q/delta'][...]
    q['gam'] = f['q/gam'][...]
    
    f.close()

    return (r, n, phi, q)
if __name__ == "__main__":
    main()
