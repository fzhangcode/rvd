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

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27
import rvd_var

##pool =None
pool = mp.Pool(processes=5)
tocfilename = "../../data/synthetic_toc_p1.txt"
toc = pd.read_table(tocfilename)

seed = 199096

# Estimate the model for the cases
logging.debug("Processing control data.")
h5FileName = "Control.hdf5"
try:
    with h5py.File(h5FileName, 'r') as f:
        pass
except IOError as e:
    controlFileList = ["../2013-08-06_Downsample_Read_Depth/depth_chart/10/%s"\
     % filename for filename in toc.Filename[toc.isRef=='Y']]
    (r, n, loc, refb) = rvd27.load_depth(controlFileList)
    (casephi, caseq) = rvd_var.ELBO_opt(r, n, seed = seed)
    logging.debug("Control data processing is done.")


# Estimate the model for the cases
for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
    logging.debug("Processing dilution: %0.1f" % dilution)
    
    h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)

    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        caseFileList = ["../2013-08-06_Downsample_Read_Depth/depth_chart/10/%s" \
        % filename for filename in toc.Filename[toc.Dilution==dilution]]
        (r, n, loc, refb) = rvd27.load_depth(caseFileList)
        logging.debug("Dilution: %0.1f processing is done." % dilution)