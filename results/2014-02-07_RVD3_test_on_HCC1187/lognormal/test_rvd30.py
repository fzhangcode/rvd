
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
T= 10/50
    
def main():
    controlFile = "./Control.hdf5"
    caseFile = "./Case.hdf5"
    outputFile='./vcf'
    rvd30.test(controlFile,caseFile,T=T,N=1000,outputFile=outputFile)
    
if __name__ == '__main__':
    main()

