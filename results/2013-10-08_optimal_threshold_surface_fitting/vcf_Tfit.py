import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

os.chdir('S:/yhe2/Research/rvd2/results/2013-10-08_optimal_threshold_surface_fitting')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    try:
        with h5py.File('./optimalT311.hdf5', 'r') as f:
            T_fit = f['T_fit'][...]        
            f.close()
            logging.debug('optimalT.hdf5 exits.')
    except IOError as e:
        logging.debug('Generate optimalT311.hdf5 first')
        sys.exit(1)
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
                  '2013-08-14_Compute_ROC_Synthetic_avg100',\
                  '2013-08-14_Compute_ROC_Synthetic_avg1000',\
                  '2013-08-14_Compute_ROC_Synthetic_avg10000')

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
    	for f in folderList:
            path='vcf/fit/311/%s' % str(10**(folderList.index(f)+1))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
            if d is not 100.0 and d is not 10.0:
                T = T_fit[dilutionList.index(d),folderList.index(f)]
            else:
                T=0               
            rvd27.test(controlFile,caseFile,T=0,N=1000,outputFile=outputFile)
if __name__ == '__main__':
    main()

