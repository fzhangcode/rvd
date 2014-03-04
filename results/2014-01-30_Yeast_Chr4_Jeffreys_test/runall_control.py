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
import pdb
# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd29 directory at front of the path

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd29

##pool =None
pool = mp.Pool(processes=20)

gibbs_nsample = 4000
mh_nsample = 5

logging.debug("Gibbs step size=%d" % gibbs_nsample)
logging.debug("mh sample size=%d" % mh_nsample)

# Estimate the model for the control
logging.debug("Processing control data.")
h5FileName = "Control.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    filename = "gen007_test.dc"
    controlFileList = ["../../../../../freeze/kvitek2013/samples_for_Chr4/DC/%s" % filename]
    (r, n, loc, refb) = rvd29.load_depth(controlFileList)
    phi, theta_s, mu_s, M_s = rvd29.mh_sample(r, n, gibbs_nsample=gibbs_nsample,mh_nsample=mh_nsample, burnin=0.2, pool=pool)
    logging.debug("Saving model in %s" % h5FileName)
    rvd29.save_model(h5FileName, phi, mu=mu_s, M=M_s, theta=theta_s, r=r, n=n, loc=loc,
           refb=refb)