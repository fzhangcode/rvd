
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
import rvd29

pool = mp.Pool(processes=40)
# T=10/Depth
T=10/400

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    controlFile = "./Control.hdf5"
    #caseFile = "./Case322.hdf5"
    #outputFile = './vcf_2'
    caseFile = "./Case133.hdf5"
    outputFile='./vcf_1'
    rvd29.test(controlFile,caseFile,T=T,N=1000,outputFile=outputFile, pool=pool)

if __name__ == '__main__':
    main()

