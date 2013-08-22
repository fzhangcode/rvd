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
import rvd27 as rvd


pool = mp.Pool(processes=48)
indexfilename = '../../data/Synthetic_BAM_files/indexFile.csv'
indexfile = pd.read_csv(indexfilename,skipinitialspace=True)


gibbs_nsample = 4000
mh_nsample = 50

logging.debug("Gibbs step size=%d" % gibbs_nsample)
logging.debug("mh sample size=%d" % mh_nsample)

# Estimate the model for the cases
logging.debug("Processing control data.")
h5FileName = "Control.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    controlFileList = ["depth_chart/%s" % filename.replace('bam','dc') for filename in indexfile.pair1Bam[indexfile.isRef=='Y']]+\
                      ["depth_chart/%s" % filename.replace('bam','dc') for filename in indexfile.pair2Bam[indexfile.isRef=='Y']]

    (r, n, loc, refb, ) =rvd.load_depth(controlFileList)
    phi, theta_s, mu_s = rvd.mh_sample(r, n, gibbs_nsample=gibbs_nsample,mh_nsample=mh_nsample, burnin=0.2, pool=pool)
    logging.debug("Saving model in %s" % h5FileName)
    rvd.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n, loc=loc,
           refb=refb)

# Estimate the model for the cases
for dilution in np.unique([float(d[8:]) for d in indexfile[indexfile.isRef=='N']['sample Name']]):
    logging.debug("Processing dilution: %0.1f" % dilution)
    
    h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)

    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        caseFileList = ["depth_chart/%s" % filename.replace('bam','dc') for filename in \
                        indexfile.pair1Bam[[float(d[8:]) for d in indexfile[indexfile.isRef=='N']['sample Name']]==dilution]]+\
                        ["depth_chart/%s" % filename.replace('bam','dc') for filename in \
                         indexfile.pair1Bam[[float(d[8:]) for d in indexfile[indexfile.isRef=='N']['sample Name']]==dilution]]
        
        (r, n, loc, refb) = rvd.load_depth(caseFileList)
        phi, theta_s, mu_s = rvd.mh_sample(r, n, gibbs_nsample=gibbs_nsample,mh_nsample=mh_nsample, burnin=0.2, pool=pool)
        logging.debug("Saving model in %s" % h5FileName)
        rvd.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n,loc=loc,
           refb=refb)
