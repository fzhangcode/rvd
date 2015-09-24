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

#rvddir = os.path.join('../../src/python/rvd27')
#sys.path.insert(0, rvddir)
import rvd2_var

##pool =None
pool = mp.Pool(processes=48)
tocfilename = "../2015-09-14_Run_rvd2_var_on_synthetic_data/synthetic_toc_p1.txt"
toc = pd.read_table(tocfilename)

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
    controlFileList = ["../2015-09-14_Run_rvd2_var_on_synthetic_data/depth_chart/10/%s" % filename for filename in toc.Filename[toc.isRef=='Y']]
    (r, n, loc, refb) = rvd2_var.load_depth(controlFileList)
    controlphi, controlq = rvd2_var.ELBO_opt(r, n, seed = 19860522, pool=48)
    logging.debug("Saving model in %s" % h5FileName)
    rvd2_var.save_model(h5FileName, r, n, controlphi, controlq, loc, refb)


# Estimate the model for the cases
for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
    logging.debug("Processing dilution: %0.1f" % dilution)
    
    h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)

    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        caseFileList = ["../2015-09-14_Run_rvd2_var_on_synthetic_data/depth_chart/10/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
        (r, n, loc, refb) = rvd2_var.load_depth(caseFileList)
        casephi, caseq = rvd2_var.ELBO_opt(r, n, seed = 19860522, pool=40)
        logging.debug("Saving model in %s" % h5FileName)
        rvd2_var.save_model(h5FileName, r, n, casephi, caseq, loc, refb) 

