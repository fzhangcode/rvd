# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
# import pandas as pd
# import pickle
# import multiprocessing as mp
import h5py
# import logging
import pdb
# import scipy.stats as ss

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

def main():
        position = ['chr7:154743899', 'chr7:154749704', 
        'chr7:154753635', 'chr7:154754371', 'chr7:154758813',
         'chr7:154760439', 'chr7:154766700', 'chr7:154766732', 
         'chr7:154777014', 'chr7:154777118', 'chr7:154780960']

	controlFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Control.hdf5'
	caseFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Case.hdf5'

        fname = 'HCC1187_partial.hdf5'

        ## Import the data

        if os.path.isfile(fname):
                with h5py.File(fname, 'r') as f:
                        muControl1 = f['muControl'][...]
                        muCase1 = f['muCase'][...]
                        muZ = f['Z'][...]
        else:
                with h5py.File(controlFile, 'r') as f:
                        muControl = f['mu'][...]
                        locControl = f['loc'][...]

                with h5py.File(caseFile, 'r') as f:
                        muCase = f['mu'][...]
                        locCase = f['loc'][...]


                idx = []

        	for pos in position:
                        if np.nonzero(locControl == pos)[0][0] is None:
                                print pos
                        else:
                		idx.append(np.nonzero(locControl == pos)[0][0])

                idx = np.array(idx)   
         
                muControl1 = muControl[idx]
                muCase1 = muCase[idx]

                N = 2000

                (muZ,_,_) =rvd27.sample_post_diff(muCase1, muControl1, N) # sample Z


                # Save the model

                h5file = h5py.File(fname, 'w')
                if muControl1 is not None:
                        h5file.create_dataset('muControl', data=muControl1,
                                chunks=True, fletcher32=True, compression='gzip')
                if muCase1 is not None:
                        h5file.create_dataset('muCase', data=muCase1,
                                chunks=True, fletcher32=True, compression='gzip')
                if muZ is not None:
                        h5file.create_dataset('Z', data=muZ,
                                chunks=True, fletcher32=True, compression='gzip')

                if position is not None:
                        h5file.create_dataset('loc', data=position,
                                chunks=True, fletcher32=True, compression='gzip')
                h5file.close()


        ## plot histogram
    
        for i in xrange(len(position)):
                fig = plt.figure(figsize=(12,3))
                ax1 = fig.add_subplot(1,2,1)
                ax1.hist(muCase1[i,:].T,20)
                ax1.hist(muControl1[i,:].T,20)
                ax1.legend( ['Case','Control'])
                ax1.set_xlabel('mu')
                plt.xticks(rotation=25) 


                ax2 = fig.add_subplot(1,2,2)
                ax2 = ax2.hist(muZ[i,:].T, 20)
                # ax2.legend(['Z',])
                fig.suptitle('position%s' %position[i].split(':')[1])
                plt.xticks(rotation=25) 
                plt.savefig('position%s.png' %position[i].split(':')[1])
                plt.tight_layout()

        ## Test if these positions will be called
        # pdb.set_trace()

        # roi = [(0.0,np.inf)]
        # p1 = rvd27.bayes_test(muZ, roi, type = 'open')

        # pdb.set_trace()

        # roi = [(-np.inf,0.0)]
        # p2 = rvd27.bayes_test(muZ, roi, type = 'open')

        # pdb.set_trace()
if __name__ == "__main__":
	main()
