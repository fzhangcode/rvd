
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss
import multiprocessing as mp

# Insert the source directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd27

pool = mp.Pool(processes=40)
# T=10/Depth
T= 0

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    
    seed = 19860522
    a = 0.05 
    
    # call the control and case from ../2014-06-10_call_mutations_in_MTH1_by_RVD3/improper/
    controlFile = "./control_HCC1187BL.hdf5"
    caseFile = "./case_HCC1187C.hdf5"
    #outputFile="./mutations_extend_PAXIP1_somatic" 
    outputFile="./mutations_extend_PAXIP1_germline" 
        
    #rvd27.test(controlFile,caseFile,T=T,N=1000,outputFile=outputFile)
    #rvd27.somatic_test(controlFile, caseFile, N=1000, intvl=[0, np.inf], alpha = a, chi2 = True, outputFile=outputFile, seedint = seed)
    rvd27.germline_test(controlFile, intvl=[0.05, np.inf], alpha = a, chi2 = True, outputFile=outputFile)


if __name__ == '__main__':
    main()

