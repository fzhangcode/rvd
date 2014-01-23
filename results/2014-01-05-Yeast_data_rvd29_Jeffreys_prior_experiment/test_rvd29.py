
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
import rvd29

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    folderList = ('./2014-01-05-Yeast_data_rvd29_Jeffreys_prior_experiment/')

        for f in folderList:
            #path='vcf/%s' % str(100**(folderList.index(f))+1)
            path='vcf/%s' % str(10000*(folderList.index(f)))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "./%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "./%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
            rvd29.test(controlFile,caseFile,T=0,N=1000,outputFile=outputFile)

if __name__ == '__main__':
    main()

