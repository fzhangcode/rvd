# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

import h5py
import pdb


# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd28')
sys.path.insert(0, rvddir)
import rvd27

def main():

        ## all positions called by muTect (for the notebook). 
        # codes in below doesn't fit these positions and 
        # needs to be taylored for the Notebook.

        ## all somatic positions called by RVD2 differnce test.
        #position = ['chr7:154749704', 'chr7:154753635', 'chr7:154754371',
        # 'chr7:154758813','chr7:154760439','chr7:154766732',
        # 'chr7:154766832','chr7:154777118']
         
        position = ['chr7:154749704', 'chr7:154753635', 'chr7:154754371',
         'chr7:154758813','chr7:154760439','chr7:154766732',
         'chr7:154766832','chr7:154777118', 'chr7:154787612', 'chr7:154788798']

	#controlFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Control.hdf5'
	#caseFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Case.hdf5'
        #controlFile='../../../../../yhe2/Research/rvd2/results/2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Control.hdf5'
        #caseFile='../../../../../yhe2/Research/rvd2/results/2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Case.hdf5'
        
        controlFile = '../2015-02-11_Expand_region_of_interest_of_PAXIP1_gene/control_HCC1187BL.hdf5'
        caseFile = '../2015-02-11_Expand_region_of_interest_of_PAXIP1_gene/case_HCC1187C.hdf5'
        
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
        drotation = 15
        legendsize = 11
              
        ## plot the histogram of the first three positions as subplots in a figure.
        J = 3
        fig = plt.figure(figsize=(12,3*J))
        for i in xrange(J):                
                ax1 = fig.add_subplot(J,2,2*i+1)
                ax1.hist(muCase1[i,:].T,20)
                ax1.hist(muControl1[i,:].T,20)
                ax1.legend( ['Case','Control'],prop={'size':legendsize})
                plt.xticks(rotation=drotation)
                ax1.set_xlabel(r'$\mu$') 
                ax1.set_ylabel('frequency') 
                ax1.text(1.1, 1.1,'Position chr7:%s' %position[i].split(':')[1], 
                    horizontalalignment='center', 
                    verticalalignment='center',
                    transform=ax1.transAxes)


                ax2 = fig.add_subplot(J,2,2*i+2)
                ax2.hist(muZ[i,:].T, 20)
                ax2.legend(['Case-Control',],prop={'size':legendsize})
                plt.xticks(rotation=drotation) 
                ax2.set_xlabel(r'$\Delta \mu$') 

        plt.tight_layout(rect=[0, 0, 0.99, 0.99])
        plt.savefig('histogram_neg.pdf')

        J =  len(position) - 3
        fig = plt.figure(figsize=(12,3*J))

        for i in xrange(J):            
                ax1 = fig.add_subplot(J,2,2*i+1)
                ax1.hist(muCase1[i+3,:].T,20)
                ax1.hist(muControl1[i+3,:].T,20)
                ax1.legend( ['Case','Control'],
                    prop={'size':legendsize})
                plt.xticks(rotation=drotation)
                ax1.set_xlabel(r'$\mu$') 
                ax1.set_ylabel('frequency') 
                ax1.text(1.1, 1.1,'Position chr7:%s' %position[i+3].split(':')[1], 
                    horizontalalignment='center', 
                    verticalalignment='center',
                    transform=ax1.transAxes)

                ax2 = fig.add_subplot(J,2,2*i+2)
                ax2.hist(muZ[i+3,:].T, 20)
                ax2.legend(['Case-Control',],prop={'size':legendsize})
                plt.xticks(rotation=drotation) 
                ax2.set_xlabel(r'$\Delta \mu$') 

        plt.tight_layout(rect=[0, 0, 0.99, 0.99])
        plt.savefig('histogram_pos.pdf')

if __name__ == "__main__":
	main()
