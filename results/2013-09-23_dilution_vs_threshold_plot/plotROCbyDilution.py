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
    dilutionList = [0.1,0.3,1.0,10.0,100.0]

    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
                  '2013-08-14_Compute_ROC_Synthetic_avg100',\
                  '2013-08-14_Compute_ROC_Synthetic_avg1000',\
                  '2013-08-14_Compute_ROC_Synthetic_avg10000')

    OPTT=[]
    COV=[]
    for f in folderList:
        cov=[]
        optT=[]
        for d in dilutionList:            
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}

            [_,_,cov_s,_,optt_s] = ROCpoints(controlFile,caseFile, d)
            cov.append(cov_s)
            optT.append(optt_s)
        COV.append(cov)
        OPTT.append(optT)
    COV=np.mean(COV,0)
    print COV
    print OPTT

    OPTT=np.array(OPTT)
    dilutionList=np.array(dilutionList)

    fig=plt.figure()
    
        

def ROCpoints(controlFile,caseFile, d,N=1000, P=0.95, chi2=False):
    # Load the model samples
    (controlPhi, controlTheta, controlMu, controlLoc, controlR, controlN) = rvd27.load_model(controlFile)
    
    (casePhi, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd27.load_model(caseFile)

    # Extract the common locations in case and control
    caseLocIdx = [i for i in xrange(len(caseLoc)) if caseLoc[i] in controlLoc]
    controlLocIdx = [i for i in xrange(len(controlLoc)) if controlLoc[i] in caseLoc]
    
    caseMu = caseMu[caseLocIdx,:]
    controlMu = controlMu[controlLocIdx,:]
    caseR = caseR[:,caseLocIdx,:]
    controlR = controlR[:,controlLocIdx,:]
    caseN = caseN[:,caseLocIdx]
    controlN = controlN[:,controlLocIdx]
    
    loc = caseLoc[caseLocIdx]
    J = len(loc)
    pos = np.arange(85,346,20)
    posidx = [i for i in xrange(J) if int(loc[i][8:]) in pos]

    
    # Sample from the posterior Z = muCase - muControl        
    (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu-casePhi['mu0'], controlMu-controlPhi['mu0'], N)
    
    # Compute cumulative posterior probability for regions (Threshold,np.inf)
    T = np.linspace(np.min(np.min(Z)), np.max(np.max(Z)), num=300)
    pList = [rvd27.bayes_test(Z, [(t, np.inf)]) for t in T]
    
    # mutation classification
    clsList = np.array((np.array(pList)>P).astype(int))
    clsList = clsList.reshape((clsList.shape[0],clsList.shape[1]))# category list

    
    # chi2 test for goodness-of-fit to a uniform distribution for non-ref bases
    if chi2:
        nRep = caseR.shape[0]
        chi2Prep = np.zeros((J,nRep))
        chi2P = np.zeros(J)
        for j in xrange(J):
                chi2Prep[j,:] = np.array([rvd27.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
                if np.any(np.isnan(chi2Prep[j,:])):
                    chi2P[j] = 1
                else:
                   chi2P[j] = 1-ss.chi2.cdf(-2*np.sum(np.log(chi2Prep[j,:] + np.finfo(float).eps)), 2*nRep) # combine p-values using Fisher's Method
        
        clsList2 = np.array((np.array(chi2P)<0.05/J).astype(int))        
        clsList2 = np.tile(clsList2,(clsList.shape[0],1))
        clsList = np.array(((clsList+clsList2)==2).astype(int))
    
    # false postive rate
    fpr = np.array([float(sum(clsList[i])-sum(clsList[i,np.array(posidx)]))/(clsList.shape[1]-len(posidx)) for i in xrange(clsList.shape[0])])
    
    # true positive rate
    tpr = np.array([float(sum(clsList[i,np.array(posidx)]))/len(posidx) for i in xrange(clsList.shape[0])])

    cov = np.median(caseN)

    # return information for mu bar plot at called positions under optimal threshold.
    distance=np.sum(np.power([fpr,tpr-1],2),0)
    Tidx=distance.argmin()

    optimalT=T[Tidx]
    return fpr,tpr, cov, T, optimalT

if __name__ == '__main__':
    main()

