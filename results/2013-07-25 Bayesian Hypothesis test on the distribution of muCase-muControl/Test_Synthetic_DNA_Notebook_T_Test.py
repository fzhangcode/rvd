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
import scipy
import pdb
print os.getcwd()
import scipy.stats as ss

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd28')
sys.path.insert(0, rvddir)
import rvd28

dilution_opt=(0.1,0.3,1.0,10.0,100.0)
gibbs_nsample_opt=[400]
##mh_nsample_opt=[10,50,100,1000]
mh_nsample_opt=[50,100]


var=False # whether the variance of muCase_s and muControl_s is assumed equal or not

for n in xrange(len(gibbs_nsample_opt)):
    for m in xrange(len(mh_nsample_opt)):
        controlFile="ngibbs="+str(gibbs_nsample_opt[n])+"_nmh="+str(mh_nsample_opt[m])+"_Control.hdf5"
        with h5py.File(controlFile, 'r') as f:
            muControl_s = f['mu'][...]
        (Nmh,ngibbs,J)=np.shape(muControl_s)

        mu_percent=[]
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

            with h5py.File(caseFile, 'r') as f:
                muCase_s = f['mu'][...]

            Z_mh=muCase_s-muControl_s
            Z_gibbs=np.median(muCase_s-muControl_s,0)

            ## plot histogram
            position1=7
            fig=plt.figure(figsize=(19,10))
##            title='ngibbs='+str(gibbs_nsample_opt[n])+'_nmh='+str(mh_nsample_opt[m])+'dilution'+str(dilution_opt[d])
##            fig.title(title)
            bins=40
            ax1=fig.add_subplot(3,2,1)
            ax1.hist(Z_gibbs[:,position1],bins)
            a='Hist(Z) when ngibbs='+str(gibbs_nsample_opt[n])+'_nmh='+str(mh_nsample_opt[m])+' position='
            title='position='+str(position1+1)
            ax1.set_title(title)

            position2=64
            ax2=fig.add_subplot(3,2,2)
            ax2.hist(Z_gibbs[:,position2],bins)
            title='position='+str(position2+1)
            ax2.set_title(title)

            position3=106
            ax3=fig.add_subplot(3,2,3)
            ax3.hist(Z_gibbs[:,position3],bins)
            title='position='+str(position3+1)
            ax3.set_title(title)


            position4=204
            ax4=fig.add_subplot(3,2,4)
            ax4.hist(Z_gibbs[:,position4],bins)
            title='position='+str(position4+1)
            ax4.set_title(title)

            position5=224
            ax5=fig.add_subplot(3,2,5)
            ax5.hist(Z_gibbs[:,position5],bins)
            title='position='+str(position5+1)
            ax5.set_title(title)

            position6=84
            ax6=fig.add_subplot(3,2,6)
            ax6.hist(Z_gibbs[:,position6],bins)
            title='position='+str(position6+1)
            ax6.set_title(title)           

            title='hist(Z) when nGibbs='+str(gibbs_nsample_opt[n])+' nMH='+str(mh_nsample_opt[m])+'dilution='+str(dilution_opt[d])+'.pdf'
 
            plt.savefig(title)

            # compute the percentage of mu that is larger than 0; Actually consider should compute the percecntage of theta higher than 0
            mu_plus=[]
            for j in xrange(J):
                indices_plus=[k for k,v in enumerate(Z_gibbs[:,j]) if v>0]
                percent_temp=np.sum(Z_gibbs[indices_plus,j])/(np.sum(np.absolute(Z_gibbs[:,j])))
                mu_plus.append(percent_temp)

            mu_plus=np.array(mu_plus)
            mu_percent.append(mu_plus)
        mu_percent=np.array(mu_percent)
        
        pos=np.arange(J)
        mpos=np.arange(85,346,20)
        
        fig2=plt.figure(figsize=(9,16))
        for d in xrange(5):
            ax=fig2.add_subplot(3,2,d+1)
            ax.plot(pos,mu_percent[d,:],marker='o')
            ax.plot(mpos-1,mu_percent[d,mpos-1],color='r',marker='o',linestyle='None')

            ax.set_xlabel('Location')
            ax.set_ylabel('p')
            ax.set_title('dilution='+str(dilution_opt[d]))
        title='mu_p when ngibbs='+str(gibbs_nsample_opt[n])+'_nmh='+str(mh_nsample_opt[m])+'.jpg'            
        plt.savefig(title)

