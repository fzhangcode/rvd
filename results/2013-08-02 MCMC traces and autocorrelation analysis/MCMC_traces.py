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

                title="Case_dilution="+str(dilution)+"_position="+str(position)+"_replicate="+str(replicate)
                title = title.replace(".", "_", 1)
                MCMCtraces(caseFile,position,replicate,dilution=dilution,title=title,figform=figform)

            controlFile="ngibbs="+str(ngibbs)+"_nmh="+str(nmh)+"_Control.hdf5"
            title="Control_position="+str(position)+"_replicate="+str(replicate)
            MCMCtraces(controlFile,position,replicate,dilution=None,title=title,figform=figform)

    
def MCMCtraces(h5FileName,position,replicate,dilution=None,title=None,figform='.pdf'):


    with h5py.File(h5FileName, 'r') as f:
        mu_s = f['mu'][...]
        theta_s = f['theta'][...]
        
    fig=plt.figure(figsize=(18,9))
    
    if dilution is not None:
        plt.suptitle('Case dilution='+str(dilution)+' position='+str(position)+' replicate='+str(replicate))
    else:
        plt.suptitle('Control position='+str(position)+' replicate='+str(replicate))
    
    '''subplot for theta'''
    (N,J,G)=np.shape(theta_s)
    ax1=fig.add_subplot(2,3,1)
    
    ax1.plot(np.arange(G),theta_s[replicate-1,position-1,:])
    ax1.set_title('MCMC trace for theta')
    ax1.set_xlabel('Gibbs steps')
    ax1.set_ylabel('theta')

    ax2=fig.add_subplot(2,3,2)
    ax2.acorr(theta_s[replicate-1,position-1,:], usevlines=True, normed=True, detrend=matplotlib.mlab.detrend_mean, maxlags=1599, lw=1)
    ax2.set_title('MCMC autocorr for theta')
    ax2.set_xlabel('lags')
    ax2.set_ylabel('autocorr')

    '''lag1 plot'''
    ax3=fig.add_subplot(2,3,3)
    ax3.plot(theta_s[replicate-1,position-1,0:G-2],theta_s[replicate-1,position-1,1:G-1],'*')
    ax3.set_title('lag1 plot of theta')
    ax3.set_xlabel('Yi-1')
    ax3.set_ylabel('Yi')
    
    '''subplot for mu'''             
    ax4=fig.add_subplot(2,3,4)
    ax4.plot(np.arange(G),mu_s[position-1,:])
    ax4.set_title('MCMC trace for mu')
    ax4.set_xlabel('Gibbs steps')
    ax4.set_ylabel('mu')
        

    ax5=fig.add_subplot(2,3,5)
    ax5.acorr(mu_s[position-1,:], usevlines=True, normed=True,detrend=matplotlib.mlab.detrend_mean, maxlags=1599, lw=2)
    ax5.set_title('MCMC autocorr for mu')
    ax5.set_xlabel('lags')
    ax5.set_ylabel('autocorr')


    ax6=fig.add_subplot(2,3,6)
    ax6.plot(mu_s[position-1,0:G-2],mu_s[position-1,1:G-1],'*')
    ax6.set_title('lag1 plot of mu')
    ax6.set_xlabel('Yi-1')
    ax6.set_ylabel('Yi')

                  
    plt.savefig(title+figform)



             
             
if __name__ == '__main__':
    main()

