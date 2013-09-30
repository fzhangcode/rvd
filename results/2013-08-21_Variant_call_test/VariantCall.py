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
    dilutionList = (0.3,)    
    folder = '2013-08-14_Compute_ROC_Synthetic_avg10'
    N=1000 # Z sampling size  

    controlFile = "../%s/Control.hdf5" %folder
  
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
        caseFile = "../%(folder)s/%(file)s"%{'folder':folder,'file':caseFile}  
        [Loc, caseMu, controlMu, postP, chi2P, call] = rvd27.test(controlFile, caseFile, T=0.00195, N=10000,outputFile=None)

    Loc = Loc.compress(call,axis=0)

    if len(Loc) is not 0:
        CallCaseMu = caseMu.compress(call,axis=0)
        CallControlMu = controlMu.compress(call,axis=0)

        ind = np.arange(CallCaseMu.shape[0])
        width=0.35

        fig, ax = plt.subplots(figsize=(11,8))
        rects1 = ax.bar(ind, np.mean(CallControlMu,1), width, color='0.8', yerr=np.std(CallControlMu,1))
        rects2 = ax.bar(ind+width, np.mean(CallCaseMu,1), width, color='0.5', yerr=np.std(CallCaseMu,1))

        # add some
        ax.set_ylabel('Minor Allele Frequency')
        ax.set_xlabel('Called Locations')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( [str(x)[8:] for x in Loc] )
        ax.set_ylim
        ax.legend( (rects1[0], rects2[0]), ('Control', 'Case') )
        
        plt.savefig('Variant_Call.png')
    else:
        print 'No variant is called'

  
if __name__ == '__main__':
    main()

