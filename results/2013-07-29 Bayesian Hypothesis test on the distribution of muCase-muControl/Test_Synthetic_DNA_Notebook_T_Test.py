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
import pyroc

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

##format='%(levelname)s:%(module)s:%(message)s')

def main():
    dilution_opt=(0.1,0.3,1.0,10.0,100.0)
    gibbs_nsample_opt=[4000]
    ##mh_nsample_opt=[10,50,100,1000]
    mh_nsample_opt=[50]
    N=1000
    steps=100
    figform='.png'
    
    
    
    for n in xrange(len(gibbs_nsample_opt)):
        for m in xrange(len(mh_nsample_opt)):
            ngibbs=gibbs_nsample_opt[n]
            nmh=mh_nsample_opt[m]
            '''ROC plot''' 
            figroc=plt.figure(figsize=(8,13))
            titleroc='ROC plot when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+' steps='+str(steps)
            plt.suptitle(titleroc)
            for d in xrange(5):
                dilution=dilution_opt[d]

                h5FileName='ROC points when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+' dilution='+str(dilution)+' steps='+str(steps)+'.hdf5'

                try:                   
                    with h5py.File(h5FileName, 'r') as f:
                        logging.debug("Z and ROC points dataset already exist")
                        muControl_s = f['muControl_s'][...]
                        muCase_s = f['muCase_s'][...]
                        Z = f['Z'][...]
                        points = f['points'][...]
                    
                except IOError as e:
                    logging.debug("Processing to get Z and ROC points")
                    logging.debug("Processing dilution: %0.1f" % dilution)
                    controlFile="ngibbs="+str(ngibbs)+"_nmh="+str(nmh)+"_Control.hdf5"
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
                    '''Bayesian hypothesis testing for Z'''
                    points=BayTest(Z,steps)
                    logging.debug("Saving model in %s" % h5FileName)
                    save_rocPoints(h5FileName,muCase_s,muControl_s,Z,points)

                    
                '''plot the histograms'''
                fig1=plt.figure(figsize=(16,9))
                bins=40
                position=[8,65,106,387,205,225]
                title='hist(Z) when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+' dilution='+str(dilution)+' steps='+str(steps)
                
                subHistPlotZ(fig1,position,Z,title,figform=figform)
                
                
                title='hist(XY) when nGibbs='+str(ngibbs)+' nMH='+str(nmh)+' dilution='+str(dilution)+' steps='+str(steps)
                fig2=plt.figure(figsize=(16,9))
                subHistPlotXY(fig2,position,muCase_s,muControl_s,title,figform=figform)
    
                
                ROCsubplot(figroc,points,d+1,'dilution='+str(dilution))
            figroc.savefig(titleroc+figform)        


def subHistPlotZ(fighandle,position,Z,title=None,bins=40,subplotsize=[3,2],figform='.pdf'):
    plt.suptitle(title)
    for i in xrange(len(position)):        
        p=position[i]
        ax=fighandle.add_subplot(subplotsize[0],subplotsize[1],i+1)
        ax.hist(Z[p-1,:],bins)
        subtitle='position='+str(p)
        ax.set_title(subtitle)
    
    title=title+figform
    plt.savefig(title)

        
def subHistPlotXY(fighandle,position,muCase_s,muControl_s,title=None,bins=40,subplotsize=[3,2],figform='.pdf'):
    plt.suptitle(title)
    for i in xrange(len(position)):
        p=position[i]
        ax=fighandle.add_subplot(subplotsize[0],subplotsize[1],i+1)
        ax.hist(muCase_s[p-1,:],bins)
        ax.hist(muControl_s[p-1,:],bins)
        subtitle='position='+str(p)
        ax.set_title(subtitle)
        ax.legend( ['Case','Control'])
    title=title+figform
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

def BayTest(Z,steps):
##    if Oneside:
##        else:
##
    (J,N)=np.shape(Z)

    '''true classification'''
    TC=np.zeros((J))
    pos=np.arange(85,346,20)
    TC[pos-1]=np.ones(len(pos))

    '''thresholds'''
    step=(np.amax(Z)-np.amin(Z))/steps
    TH=np.arange(np.amin(Z),np.amax(Z)+step,step) 
    THlen=len(TH)

    '''points in ROC plot'''
    points=np.zeros((THlen,4))# [fpr,tpr,TH,p]
    
    for i in xrange(THlen):
        C=np.zeros((J))
        if i%10==0:
            print i
        th=TH[i]
        for j in xrange(J):
            indices=[k for k,v in enumerate(Z[j,:]) if v>th]

            p=np.sum(np.absolute(Z[j,indices]))/(np.sum(np.absolute(Z[j,:])))
            if p==.5:
                if np.random.choice([True, False]):
                    C[j]=1
            else:
                C[j]=int(p>.5)

        fpC=np.copy(C)
        fpC[pos-1]=np.zeros(len(pos))
        fpr=float(sum(fpC))/(len(TC)-len(pos))

        tpr=float(sum(C[pos-1]))/len(pos)

        '''accuracy'''
        acc=float(sum(map(int,C==TC)))/len(TC)
        points[i,:]=[fpr,tpr,th,acc]
        
    return points

    '''ROC plot'''    
def ROCsubplot(figroc,points,NOsubplot,subtitle):
    axx=figroc.add_subplot(3,2,NOsubplot)
    axx.plot(points[:,0],points[:,1],color='b',marker='o')

    axx.set_title(subtitle)
    axx.grid(True)
    axx.set_ylabel('True Positive Rate')
    axx.set_xlabel('False Positive Rate')

    '''text the threshold when the ROC point is closest to the (0,1)point'''    
    dpoints=[points[:,0],points[:,1]-1]
    dpoints=np.power(dpoints,2)
    distance=np.sum(dpoints,0)
    ind=distance.argmin()
    axx.plot(points[ind,0],points[ind,1],color='r',marker='o')
    axx.text(0.3,0.6,'TH='+str(points[ind,2]))
    axx.text(0.3,0.5,'ACC='+str(points[ind,3]))

def save_rocPoints(h5FileName,muCase_s=None,muControl_s=None,Z=None,points=None):

    # TODO add attributes to hdf5 file
    h5file = h5py.File(h5FileName, 'w')
    
    # Save the latent variables if available.
    if muCase_s is not None:
        h5file.create_dataset('muCase_s', data=muCase_s, 
                              chunks=True, fletcher32=True, compression='gzip')
    if muControl_s is not None:
        h5file.create_dataset('muControl_s', data=muControl_s, 
                              chunks=True, fletcher32=True, compression='gzip')
    
    if Z is not None:
        h5file.create_dataset('Z', data=Z, 
                              chunks=True, fletcher32=True, compression='gzip')
    if points is not None:
        h5file.create_dataset('points', data=points, 
                              chunks=True, fletcher32=True, compression='gzip')

    h5file.close()



if __name__ == '__main__':
    main()

