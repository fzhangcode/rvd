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

    folderList = ('one_replicate_synthetic_avg10',\
                  'one_replicate_synthetic_avg100',\
                  'one_replicate_synthetic_avg1000',\
                  'one_replicate_synthetic_avg10000')

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
    	for f in folderList:
            path='onerep_vcf_mu0_chi2_05/%s' % str(10**(folderList.index(f)+1))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s.vcf' %{'path':path,'dilution':str(d).replace('.','_')}
           
            rvd27.somatic_test(controlHDF5Name=controlFile,caseHDF5Name=caseFile, tau = 0,N=1000,outputFile=outputFile)
if __name__ == '__main__':
    main()

