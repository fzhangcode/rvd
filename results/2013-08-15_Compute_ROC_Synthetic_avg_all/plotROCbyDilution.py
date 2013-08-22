import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
              '2013-08-14_Compute_ROC_Synthetic_avg100',\
              '2013-08-14_Compute_ROC_Synthetic_avg1000', \
	      '2013-08-14_Compute_ROC_Synthetic_avg10000')

    
    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)

    i=0
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        i=i+1
        plt.figure(i)
        label=[]
        for f in folderList:            
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
             # ROC
            [fpr,tpr,cov] = ROCpoints(controlFile,caseFile)            
            
            plt.plot(fpr,tpr, label='%0.1f' % cov)

        plt.legend(loc=4)            
        plt.title('%0.1f%% Mutant Mixture' % d)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

        title = 'ROC_dilution%0.1f' % d
        plt.savefig('%s.pdf' % title.replace('.','_',1))

def ROCpoints(controlFile,caseFile,N=1000,P=0.95):
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
    
    loc = caseLoc[caseLocIdx]
    J = len(loc)
    #pos = np.arange(85,346,20)
    idx = (384, 7, 29, 51, 73, 95, 118, 140, 162, 184, 206, 229, 251, 273)
    #pos = [loc[i] for i in xrange(J) if loc[i] in pos]

    # Sample from the posterior Z = muCase - muControl        
    (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu-casePhi['mu0'], controlMu-controlPhi['mu0'], N)
    
    # Compute cumulative posterior probability for regions (Threshold,np.inf) 
    pList = [rvd27.bayes_test(Z, [(T, np.inf)]) for T in np.linspace(np.min(np.min(Z)), np.max(np.max(Z)), num=300)]

    # mutation classification
    clsList = np.array((np.array(pList)>P).astype(int))
    clsList = clsList.reshape((clsList.shape[0],clsList.shape[1]))# category list
    
    # false postive rate
    fpr = np.array([float(sum(clsList[i])-sum(clsList[i,idx]))/(clsList.shape[1]-len(idx)) \
           for i in xrange(clsList.shape[0])])
    
    # true positive rate
    tpr = np.array([float(sum(clsList[i, idx]))/len(idx) for i in xrange(clsList.shape[0])])
   
    cov = np.median(caseN)
    return fpr,tpr, cov

if __name__ == '__main__':
    main()

