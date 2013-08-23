import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import argparse

##os.chdir('S://yhe2//Research//rvd2//results//2013-08-21_depth_M_plot')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
##    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    dilutionList = (1.0,)
    
    folder = '2013-08-19_Compute_ROC_Synthetic_avg10000'
    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)
    controlFile = "../%s/Control.hdf5" %folder
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
        caseFile = "../%(folder)s/%(file)s" %{'folder':folder,'file':caseFile}
        T=0.0005
        outputFile='test%s.hdf5'% str(d).replace(".","_")
##        rvd27.test([controlFile,caseFile,N,T,outputFile])
        args=(controlFile, caseFile, N, T, outputFile)
        rvd27.test(args)
if __name__ == '__main__':
    main()

