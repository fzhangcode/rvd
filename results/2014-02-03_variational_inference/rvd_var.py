# -*- coding: utf-8 -*-
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



def main():
    log_level = logging.DEBUG # default logging level
    logging.basicConfig(level=log_level,
        format='%(levelname)s:%(module)s:%(message)s')
    
    (y, x, theta) = gen1()
    
    pool = None
    #pool = mp.Pool(processes=2)
    
    (phi, q) = varlap_opt(y, K=2, seed=10241978, pool=pool)

    save_model('model.hdf5', y, phi, q)
    

    

## compute sufficient statistics
def EqlogTheta(alpha, beta):
	return psi(alpha) - psi(alpha + beta)

def Eqlog1_Theta(alpha, beta):
	return psi(beta) - psi(alpha + beta)

def EqMu(delta, gamma):
	return delta / (delta + gamma) # eps?

def EqlogMu(delta, gamma):
	return psi(delta) - psi(delta + gamma)

def Eqlog1_Mu(delta, gamma):
	return psi(gamma) - psi(delta + gamma)

def EqlogGamma(delta, gamma, M):
	# Eqlog\Gamma(muM)
	N = 7
	mu = range(0, 1 , 2^N+1)
	y = [kernel(delta, gamma, x, M) for x in mu]

	return integrate.romb(y, dx = 0.1/(2^N+1))
	

def kernel(delta, gamma, mu, M):
		return ss.beta.pdf(mu, delta, gamma)*betaln(mu*M, (1-mu)*M)


## compute entropy
def BetaEntropy(a, b):
	# To compute EqlogQmu and EqlogQtheta
	#return betaln(delta,gamma) - (delta-1)psi(delta) - (gamma - 1) psi(gamma) +(delta + gamma -2) psi(delta + gamma)
	return betaln(a,b) - (a-1)psi(a) - (b - 1) psi(b) +(a + b -2) psi(a + b)


## compute ELBO
def ELBO(r, n, theta, mu, M, mu0, M0, alpha, beta, delta, gamma):
	(J,N) = n.shape

	Mu = [EqMu(delta[j], gamma[j]) for j in xrange(J)]
	logMu = [EqlogMu(delta[j], gamma[j]) for j in xrange(J)]
	log1_Mu = [Eqlog1_Mu(delta[j], gamma[j]) for j in xrange(J)]
	logTheta = [EqlogTheta(alpha[j,i], beta[j,i]) for j in xrange(J) for i in xrange(N)]
	log1_Theta = [Eqlog1_Theta(alpha[j,i], beta[j,i]) for j in xrange(J) for i in xrange(N)]

	EqlogPr = 0.0
	for j in xrange(J):
		for i in xrange(N):
			EqlogPr += - betaln(r[j,i] + 1, n[j,i] - r[j,i] +1)
			EqlogPr += r[j,i]*logTheta[j,i] + (n[j,i] - r[j,i]) * logTheta[j,i]

	EqlogPtheta = 0.0
	for j in xrange(J):
		EqlogPtheta += - betaln(mu[j]*M[j] + np.finfo(float).eps, (1-mu[j]*M[j]) + np.finfo(float).eps)
		for i in xrange(N):
			EqlogPtheta += (M[j]* Mu[j]- 1)*logTheta[j,i] +\
			(M[j]*(1 - Mu[j]) - 1)*log1_Theta[j,i]

	EqlogPmu = -betaln(mu0*M0, (1-mu0)*M0)
	for j in xrange(J):
		EqlogPmu += (M0*mu0-1)*logMu[j] + (M0*(1-mu0)-1)*log1_Mu[j]

	EqlogQtheta = 0.0
	for j in xrange(J):
		for i in xrange(N):
			EqlogQtheta += BetaEntropy(alpha[j,i], beta[j,i])

	EqlogQmu = 0.0
	for j in xrange(J):
			EqlogQmu += BetaEntropy(delta[j], gamma[j])




if __name__ == "__main__":
    main()
