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

def Eqlog1_mu(delta, gamma):
	return psi(gamma) - psi(delta + gamma)

def EqlogGamma(delta, gamma, M):
	# Eqlog\Gamma(muM)
	N = 6
	mu = range(0, 1 , 2^N+1)
	y = [kernel(delta, gamma, x, M) for x in mu]

	return integrate.romb(y, dx = 0.1/(2^N+1))
	

def kernel(delta, gamma, mu, M):
		return ss.beta.pdf(mu, delta, gamma)*betaln(mu*M, (1-mu)*M)

## compute entropy
def BetaEntropy(a, b):
	#return betaln(delta,gamma) - (delta-1)psi(delta) - (gamma - 1) psi(gamma) +(delta + gamma -2) psi(delta + gamma)
	return betaln(a,b) - (a-1)psi(a) - (b - 1) psi(b) +(a + b -2) psi(a + b)

## compute ELBO
def ELBO():
	EqlogPr = 0
	EqlogPtheta = 
	EqlogPmu = 
	EqlogQtheta = 
	EqlogQmu = 




if __name__ == "__main__":
    main()
