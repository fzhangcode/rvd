
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the source directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd30

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

    
def main():
    controlFile = "./1/Control.hdf5"
    caseFile = "./1/Case0.hdf5"
    outputFile='./vcf'
    rvd30.test(controlFile,caseFile,T=0,N=1000,outputFile=outputFile)
    
if __name__ == '__main__':
    main()

