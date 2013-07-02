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
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

# <codecell>

def load_depth(dcFileNameList):
    r=[]; n=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    for dcFileName in dcFileNameList:
        with open(dcFileName, 'r') as dcFile:
            header = dcFile.readline().strip()
            dc = dcFile.readlines()
            dc = [x.strip().split("\t") for x in dc]
    
            loc = map(int, [x[2] for x in dc])
            refb = [x[4] for x in dc]
            c = [map(int, x[5:9]) for x in dc]
            c = np.array(c)
    
            (J, K) = np.shape(c)
            n1 = np.sum(c, 1)
            r1 = np.zeros(J)
            for j in xrange(0, J):
                r1[j] = n1[j] - c[j, acgt[refb[j]]]
        r.append(r1)
        n.append(n1)
    r = np.array(r)
    n = np.array(n)
    return (r,n)

# <codecell>

pool = mp.Pool(processes=60)
tocfilename = "synthetic_toc.txt"
toc = pd.read_table(tocfilename)

# <codecell>

# Estimate the model for the cases
logging.debug("Processing control data.")
h5FileName = "control.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    controlFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.isRef=='Y']]
    (r, n) = load_depth(controlFileList)
    phi, theta_s, mu_s = rvd27.mh_sample(r, n, nsample=400, burnin=0.2, pool=pool)
    logging.debug("Saving model in %s" % h5FileName)
    rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n)

# <codecell>

# Estimate the model for the cases
for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
    logging.debug("Processing dilution: %0.1f" % dilution)
    
    h5FileName = "case%0.1f.hdf5" % dilution
    h5FileName = h5FileName.replace(".", "_", 1)
    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        caseFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
        (r, n) = load_depth(caseFileList)
        phi, theta_s, mu_s = rvd27.mh_sample(r, n, nsample=500, burnin=0.2, pool=pool)
        logging.debug("Saving model in %s" % h5FileName)
        rvd27.save_model(h5FileName, phi, mu=mu_s, theta=theta_s, r=r, n=n)

