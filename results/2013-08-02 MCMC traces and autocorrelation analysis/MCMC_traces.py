import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import pickle
import multiprocessing as mp
import h5py
import logging
import pdb
import scipy.stats as ss

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilution_opt=(0.1,0.3,1.0,10.0,100.0)
    gibbs_nsample_opt=[4000]
    mh_nsample_opt=[50]
    position=225
    replicate=2
    figform='.png'    
    
    for n in xrange(len(gibbs_nsample_opt)):
        for m in xrange(len(mh_nsample_opt)):
            ngibbs=gibbs_nsample_opt[n]
            nmh=mh_nsample_opt[m]
           
            '''ROC plot''' 
            for d in xrange(5):
                dilution=dilution_opt[d]
                
                caseFile="ngibbs="+str(gibbs_nsample_opt[n])+"_nmh="+str(mh_nsample_opt[m])+'_'
                if dilution==0.1:
                    caseFile+='Case0_1.hdf5'
                elif dilution==0.3:
                    caseFile+='Case0_3.hdf5'
                elif dilution==1.0:
                    caseFile+='Case1_0.hdf5'
                elif dilution==10.0:
                    caseFile+='Case10_0.hdf5'
                else:
                    caseFile+='Case100_0.hdf5'

                title=" dilution="+str(dilution)+" position="+str(position)+" replicate="+str(replicate)+" Case"
                MCMCtraces(caseFile,position,replicate,title,figform='.pdf')

            controlFile="ngibbs="+str(ngibbs)+"_nmh="+str(nmh)+"_Control.hdf5"
            title="position="+str(position)+" replicate="+str(replicate)+" Control"
            MCMCtraces(controlFile,position,replicate,title,figform='.pdf')

    
def MCMCtraces(h5FileName,position,replicate,title=None,figform='.pdf'):
    
    def autocorr(x):
        result = np.correlate(x, x, mode='full')
        return result[result.size/2:]


    with h5py.File(h5FileName, 'r') as f:
        mu_s = f['mu'][...]
        theta_s = f['theta'][...]
        
    fig=plt.figure(figsize=(16,9))
    plt.suptitle('position='+str(position)+' replicate='+str(replicate))
    
    '''subplot for theta'''
    (N,J,G)=np.shape(theta_s)
    ax1=fig.add_subplot(2,2,1)
    ax1.plot(np.arange(G),theta_s[replicate-1,position-1,:])
    ax1.set_title('MCMC trace for theta')
    ax1.set_xlabel('Gibbs steps')
    ax1.set_ylabel('theta')

    ax2=fig.add_subplot(2,2,2)
    theta_ac=autocorr(theta_s[replicate-1,position-1,:])         
    ax2.plot(np.arange(G),theta_ac)
    ax2.set_title('MCMC autocorr for theta')
    ax2.set_xlabel('Gibbs steps')
    ax2.set_ylabel('autocorr')

    
    '''subplot for mu'''             
    ax3=fig.add_subplot(2,2,3)
    ax3.plot(np.arange(G),mu_s[position-1,:])
    ax3.set_title('MCMC trace for mu')
    ax3.set_xlabel('Gibbs steps')
    ax3.set_ylabel('mu')
        

    ax4=fig.add_subplot(2,2,4)
    mu_ac=autocorr(mu_s[position-1,:])         
    ax4.plot(np.arange(G),mu_ac)
    ax4.set_title('MCMC autocorr for mu')
    ax4.set_xlabel('Gibbs steps')
    ax4.set_ylabel('autocorr')
    plt.savefig(title+figform)
             


if __name__ == '__main__':
    main()

