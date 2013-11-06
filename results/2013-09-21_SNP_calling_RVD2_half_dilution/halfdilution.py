import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
                  '2013-08-14_Compute_ROC_Synthetic_avg100',\
                  '2013-08-14_Compute_ROC_Synthetic_avg1000',\
                  '2013-08-14_Compute_ROC_Synthetic_avg10000')

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
    	for f in folderList:
            path='vcf/%s' % str(10**(folderList.index(f)+1))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
           
            rvd27.test(controlFile,caseFile,T=d/2,N=1000,outputFile=outputFile)
            
if __name__ == '__main__':
    main()

