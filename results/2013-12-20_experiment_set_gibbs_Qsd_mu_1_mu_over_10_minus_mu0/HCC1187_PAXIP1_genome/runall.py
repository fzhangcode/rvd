# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp
import logging
import random

# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27


pool = mp.Pool(processes=5)
#pool = None

gibbs_nsample = 4000
mh_nsample = 5

logging.debug("Gibbs step size=%d" % gibbs_nsample)
logging.debug("mh sample size=%d" % mh_nsample)

def fit_model(fileNameList, h5FileName):
    fileList = ["../../2013-11-07_HCC1187_PAXIP1_genome_fa/depth_chart/%s" % filename for filename in fileNameList]
    (r, n, loc, refb) = rvd27.load_depth(fileList)
    phi, theta_s, mu_s = rvd27.mh_sample(r, n, gibbs_nsample=gibbs_nsample,mh_nsample=mh_nsample, burnin=0.2, pool=pool)
    logging.debug("Saving model in %s" % h5FileName)
    rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,
       refb=refb)

def main():
    random.seed(199096)
    # Estimate the model for the cases
    logging.debug("Processing control data.")
    controlList = ['HCC1187BL_S1.dc',]
    fit_model(controlList, 'Control.hdf5')


##    # Estimate the model for the cases
##    logging.debug("Processing case data.")
##    caseList = ['HCC1187C_S1.dc',]
##    fit_model(caseList, 'Case.hdf5')

if __name__ == "__main__":
	main()

