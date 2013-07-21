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

# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27


# <codecell>

pool = mp.Pool(processes=30)
tocfilename = "synthetic_toc.txt"
toc = pd.read_table(tocfilename)

# <codecell>
nsample_opt=[4000,40000]
for i in xrange(len(nsample_opt)):
    nsample=nsample_opt[i]
    logging.debug("SamplingSize="+str(nsample))
    # Estimate the model for the cases
    logging.debug("Processing control data.")
    h5FileName = "SamplingSize="+str(nsample)+"_Control.hdf5"
    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        controlFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.isRef=='Y']]
        (r, n, loc, refb, ee) = rvd27.load_depth(controlFileList)
        phi, theta_s, mu_s = rvd27.mh_sample(r, n, nsample=nsample, burnin=0.2, pool=pool)
        logging.debug("Saving model in %s" % h5FileName)
        rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,
               refb=refb, ee=ee)

    # <codecell>

    # Estimate the model for the cases
    for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
        logging.debug("Processing dilution: %0.1f" % dilution)
        
        h5FileName = "SamplingSize="+str(nsample)+"_Case%0.1f.hdf5" % dilution
        h5FileName = h5FileName.replace(".", "_", 1)
        try:
            with h5py.File(h5FileName, 'r') as f:
                pass
        except IOError as e:
            caseFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
            (r, n, loc, refb, ee) = rvd27.load_depth(caseFileList)
            phi, theta_s, mu_s = rvd27.mh_sample(r, n, nsample=nsample, burnin=0.2, pool=pool)
            logging.debug("Saving model in %s" % h5FileName)
            rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n,loc=loc,
               refb=refb, ee=ee)

