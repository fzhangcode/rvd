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
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('2013-11-15_six_replicates_synthetic_avg10',\
                  '2013-11-15_six_replicates_synthetic_avg100',\
                  '2013-11-15_six_replicates_synthetic_avg1000',\
                  '2013-11-15_six_replicates_synthetic_avg10000')
    outputFile='optimalT.hdf5'
    try:
        with h5py.File(outputFile, 'r') as f:
                cov = f['cov'][...]
                T_ED = f['T_ED'][...]
                T_L1 = f['T_L1'][...]
                T_MCC = f['T_MCC'][...]           
                f.close()
                logging.debug('optimalT.hdf5 exits already')
    except IOError as e:
        logging.debug('Generate optimal dataset...')
        N=1000 # Z sampling size  
        chi2=False
        times=8
        T_ED=np.zeros((times,len(dilutionList),len(folderList)))
        T_L1=np.zeros((times,len(dilutionList),len(folderList)))
        T_MCC=np.zeros((times,len(dilutionList),len(folderList)))
        cov=np.zeros((times,len(dilutionList),len(folderList)))
        for i in xrange(times):
            logging.debug('Runing the %(i)d times out of %(total)d...' % {'i':i+1, 'total':times})
            for d in dilutionList:
                logging.debug("Processing dilution: %0.1f%%" % d)
                for f in folderList:
                    path='%s' % str(10**(folderList.index(f)+1))
                    controlFile = "../%s/Control.hdf5" %f
                    caseFile = "Case%s.hdf5" % str(d).replace(".","_")
                    caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
                    m=dilutionList.index(d)
                    n=folderList.index(f)
                    (T_ED[i][m][n],T_L1[i][m][n],T_MCC[i][m][n],cov[i][m][n])=ROCpoints(controlFile,caseFile, path, d, P=0.95,chi2=chi2)
        cov=np.mean(np.median(cov,0),0)
        with h5py.File(outputFile, 'w') as f:
            f.create_dataset('cov', data=cov)
            f.create_dataset('dilutionList', data=dilutionList)
            f.create_dataset('T_ED', data=T_ED)
            f.create_dataset('T_L1', data=T_L1)
            f.create_dataset('T_MCC', data=T_MCC)
            f.close()
    
def ROCpoints(controlFile,caseFile, path, d,N=1000, P=0.95, chi2=False):
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
    posidx = np.array(posidx)
    # reference list
    refList=np.zeros(J)
    refList[posidx]=np.ones(len(posidx))

    
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

    ## compute the statistics measures
    fpr=np.zeros(clsList.shape[0])
    tpr=np.zeros(clsList.shape[0])
    tnr=np.zeros(clsList.shape[0])
    mcc=np.zeros(clsList.shape[0])

    #refList Postive
    P1 = sum(refList)
    #refList Negative
    N1 = len(refList) - P1

    for j in xrange(clsList.shape[0]):       
        #True Positive
        TP = len([i for i in range(len(refList)) if refList[i]==1 and clsList[j,i]==1])
        #True Negative
        TN = len([i for i in range(len(refList)) if refList[i]==0 and clsList[j,i]==0])
        #False Positive
        FP = len([i for i in range(len(refList)) if refList[i]==0 and clsList[j,i]==1])
        #False Negative
        FN = len([i for i in range(len(refList)) if refList[i]==1 and clsList[j,i]==0])
      
        #clsList Postive
        P2 = sum(clsList[j])
        #clsList Negative
        N2 = len(clsList[j]) - P2
        
        #Sensitivity (TPR,true positive rate)
        if TP+FN != 0:
            tpr[j] = float(TP)/(TP+FN)
        else:
            tpr[j] = np.nan
            
        #Specificity (TNR, true negative rate)
        if FP+TN != 0:
            tnr[j] = float(TN)/(FP+TN)
        else:
            tnr[j] = np.nan
            
        #FPR (1-Specificity, false negative rate)
        if FP+TN != 0:
            fpr[j] = float(FP)/(FP+TN)
        else:
            fpr[j] = np.nan

        #MCC (Matthews correlation coefficient)
        if P1*N1*P2*N2 != 0:
            mcc[j] = float(TP*TN-FP*FN)/np.sqrt(P1*N1*P2*N2)
        else:
            mcc[j]= np.nan

    cov = np.median(caseN)

    # Euclidean distance
    distance1=np.sum(np.power([fpr,tpr-1],2),0)
    Tidx1=distance1.argmin()
    
    # L1 distance
    distance2 = tpr+tnr
    Tidx2=distance2.argmax()
    
    # mcc; return the index of maximum non_NAN number
    def removenan(x):
        if np.isnan(x):
            return 0
        else:
            return x
    mcc_reform = np.array(map(removenan,mcc))
    Tidx3=mcc_reform.argmax()

    print 'ED:%d' %Tidx1
    print 'L1:%d' %Tidx2
    print 'MCC:%d' %Tidx3
    
    # return information for mu bar plot at called positions under optimal threshold.
    EDpath='vcf/ED/%s' %path
    L1path='vcf/L1/%s' %path
    MCCpath='vcf/MCC/%s' %path

    save_result(EDpath,Tidx1,d,controlFile,controlLocIdx,clsList,controlMu,caseMu,caseR,loc,J)
    save_result(L1path,Tidx2,d,controlFile,controlLocIdx,clsList,controlMu,caseMu,caseR,loc,J)
    save_result(MCCpath,Tidx3,d,controlFile,controlLocIdx,clsList,controlMu,caseMu,caseR,loc,J)

    return T[Tidx1],T[Tidx2],T[Tidx3],cov
    

def save_result(path,Tidx,dilution,controlFile,controlLocIdx,clsList,controlMu,caseMu,caseR,loc,J):
    if not os.path.exists(path):
        os.makedirs(path)
    outputFile=os.path.join(path,'vcf%s.vcf' %str(dilution).replace('.','_'))

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
        
    rvd27.write_vcf(outputFile, loc, call, refb, altb, np.mean(controlMu, axis=1),np.mean(caseMu, axis=1))

if __name__ == '__main__':
    main()
