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
import scipy.stats as ss

def main():
    dilution_opt=(0.1,0.3,1.0,10.0,100.0)
    gibbs_nsample_opt=[4000]
    ##mh_nsample_opt=[10,50,100,1000]
    mh_nsample_opt=[50]
    N=1000
    for n in xrange(len(gibbs_nsample_opt)):
        for m in xrange(len(mh_nsample_opt)):
            ngibbs=gibbs_nsample_opt[n]
            nmh=mh_nsample_opt[m]
            controlFile="ngibbs="+str(ngibbs)+"_nmh="+str(nmh)+"_Control.hdf5"
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
                    
                (muCase_s,muControl_s,Z)=Zsampling(caseFile,controlFile,N)
                
                # plot the histograms
                fig1=plt.figure(figsize=(16,9))
                bins=40
                position=[8,65,106,387,205,225]
                
                subHistPlotZ(fighandle=fig1,position=position,ngibbs=ngibbs,nmh=nmh,dilution=dilution,Z=Z)
             
                fig2=plt.figure(figsize=(16,9))
                subHistPlotXY(fighandle=fig2,position=position,ngibbs=ngibbs,nmh=nmh,muCase_s=muCase_s,dilution=dilution,muControl_s=muControl_s)

def subHistPlotZ(fighandle,position,Z,ngibbs,nmh,dilution,bins=40,subplotsize=[3,2]):
    for i in xrange(len(position)):
        
        p=position[i]
        ax=fighandle.add_subplot(subplotsize[0],subplotsize[1],i+1)
        ax.hist(Z[p-1,:],bins)
        title='Z distribution at position='+str(p)
        ax.set_title(title)
    title='hist(Z) when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+'dilution='+str(dilution)+'.pdf'    
    plt.savefig(title)

        
def subHistPlotXY(fighandle,position,muCase_s,muControl_s,ngibbs,nmh,dilution,bins=40,subplotsize=[3,2]):
    for i in xrange(len(position)):
        p=position[i]
        ax=fighandle.add_subplot(subplotsize[0],subplotsize[1],i+1)
        ax.hist(muCase_s[p-1,:],bins)
        ax.hist(muControl_s[p-1,:],bins)
        title='X Y distribution at position='+str(p)
        ax.set_title(title)
        ax.legend( ['Case','Control'])
    title='hist(XY) when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+'dilution='+str(dilution)+'.pdf'    
    plt.savefig(title)

  
def Zsampling(caseFile,controlFile,N):
    def Xsampling(File,N):
        with h5py.File(File, 'r') as f:
            mu_s = f['mu'][...]
        (J,Ngibbs)= np.shape(mu_s)
        mu_sampled=np.zeros((J,N))
        for j in xrange(J):
            for i in xrange(N):
                mu_sampled[j,i]=np.random.choice(mu_s[j,:])
        return mu_s, mu_sampled
    
    (muCase_s,muCase_sampled)=Xsampling(caseFile,N)
    (muControl_s,muControl_sampled)=Xsampling(controlFile,N)

    Z=muCase_sampled-muControl_sampled

    return muCase_s, muControl_s, Z

if __name__ == '__main__':
    main()
