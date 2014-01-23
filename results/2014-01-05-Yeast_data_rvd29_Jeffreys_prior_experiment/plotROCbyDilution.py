import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the rvd29 directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd29

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0)

    folderList = ('2013-11-18_test_rvd29_synthetic_data_dc=10',\
                    '2013-11-18_test_rvd29_synthetic_data_dc=100',\
                    '2013-11-19_test_rvd29_synthetic_data_dc=1000',\
                    '2013-11-19_test_rvd29_synthetic_data_dc=10000')
    
    N=1000 # Z sampling size  
    fig=plt.figure(figsize=(10, 9))
    chi2=False
    lstyle=('-','--','-.',':')
    lcolor=('c','r','g','b')
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        ax = fig.add_subplot(2,2,dilutionList.index(d)+1)
        label=[]
        for f in folderList:
            path='vcf/%s' % str(10**(folderList.index(f)+1))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
             # ROC
            [fpr,tpr,cov, T] = ROCpoints(controlFile,caseFile, path, d, P=0.95,chi2=chi2)
##            ax.plot(fpr,tpr, color=lcolor[folderList.index(f)], label='%d' % cov)
            ax.plot(fpr,tpr,linestyle=lstyle[folderList.index(f)], color=lcolor[folderList.index(f)], linewidth=2.0, label='%d' % cov)
            ax.plot(fpr[0],tpr[0],marker='o',markerfacecolor=lcolor[folderList.index(f)])
        l = ax.legend(loc=4,prop={'size':9},title='Read depth')
        l.get_title().set_fontsize(9)
        ax.plot([0,1],[0,1],color='k',linestyle='dashed')
        ax.set_title('%0.1f%% Mutant Mixture' % d, fontsize=10)
        ax.set_xlim((-0.03,1.03))
        ax.set_ylim((-0.03,1.03))
                 
        ax.set_xlabel('1-Specificity (FPR)',fontsize=10)
        ax.set_ylabel('Sensitivity (TPR)',fontsize=10)
        
    if chi2:
        title='ROC_with_chi2'
    else:
        title='ROC_without_chi2'
    figformat='.pdf'
    plt.savefig(title+figformat)

def ROCpoints(controlFile,caseFile, path, d,N=1000, P=0.95, chi2=False):
    # Load the model samples
    (controlPhi, controlM, controlTheta, controlMu, controlLoc, controlR, controlN) = rvd29.load_model(controlFile)
    (casePhi, caseM, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd29.load_model(caseFile)
    
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
    (Z, caseMuS, controlMuS) = rvd29.sample_post_diff(caseMu-casePhi['mu0'], controlMu-controlPhi['mu0'], N)
    
    # Compute cumulative posterior probability for regions (Threshold,np.inf)
    T = np.linspace(np.min(np.min(Z)), np.max(np.max(Z)), num=300)
    pList = [rvd29.bayes_test(Z, [(t, np.inf)]) for t in T]
    
    # mutation classification
    clsList = np.array((np.array(pList)>P).astype(int))
    clsList = clsList.reshape((clsList.shape[0],clsList.shape[1]))# category list

    
    # chi2 test for goodness-of-fit to a uniform distribution for non-ref bases
    if chi2:
        nRep = caseR.shape[0]
        chi2Prep = np.zeros((J,nRep))
        chi2P = np.zeros(J)
        for j in xrange(J):
                chi2Prep[j,:] = np.array([rvd29.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
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

    # using EL distance.
##    distance=np.sum(np.power([fpr,tpr-1],2),0) 
##    Tidx=distance.argmin()
##    print Tidx


     # Using L1 distance 
    distance = 1+tpr-fpr
    Tidx=distance.argmax()
    

    outputFile=os.path.join(path,'vcf%s.vcf' %str(d).replace('.','_'))
    
    with h5py.File(controlFile, 'r') as f:
        refb = f['/refb'][...]
        f.close()
    refb = refb[controlLocIdx]
    
    altb = []
    call=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in xrange(J):
        r = np.squeeze(caseR[:,i,:]) # replicates x bases
        
        # Make a list of the alternate bases for each replicate
        acgt_r = ['A','C','G','T']
        del acgt_r[ acgt[refb[i]] ]

        altb_r = [acgt_r[x] for x in np.argmax(r, axis=1)]

        if clsList[Tidx,i]==1:
            call.append(True)
            altb.append(altb_r[0])
        else:
            altb.append(None)
            call.append(False)
            
    rvd29.write_vcf(outputFile, loc, call, refb, altb, np.mean(caseMu, axis=1), np.mean(controlMu, axis=1))
    return fpr,tpr, cov, T

if __name__ == '__main__':
    main()

