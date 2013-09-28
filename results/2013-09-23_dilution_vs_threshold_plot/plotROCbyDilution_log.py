import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss
import math
from scipy import interpolate
from mpl_toolkits.mplot3d.axes3d import Axes3D

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    
    dilutionListopt = [[0.1,0.3,1.0,10.0,100.0],[0.1,0.3,1.0,10.0]]
    
    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
                '2013-08-14_Compute_ROC_Synthetic_avg100',\
                '2013-08-14_Compute_ROC_Synthetic_avg1000',\
                '2013-08-14_Compute_ROC_Synthetic_avg10000')
    outputfileopt = ('output5.hdf5','output4.hdf5')
    figtitleopt = ('OptimalT_vs_Dilution5_log','OptimalT_vs_Dilution4_log')
    
    plt.close('all')
    for m in xrange(2):
        dilutionList=dilutionListopt[m]
        outputFile=outputfileopt[m]
        
        try:
            with h5py.File(outputFile, 'r') as f:
                cov = f['cov'][...]
                optT = f['optT'][...]
                f.close()
        except IOError as e:
            Times=15
            optT=np.zeros((Times, 4, len(dilutionList)))
            cov=np.zeros((Times, 4, len(dilutionList)))
            for i in xrange(Times):
                for f in folderList:
                    for d in dilutionList:            
                        controlFile = "../%s/Control.hdf5" %f
                        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
                        caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
                        [_,_,cov_s,_,optt_s] = ROCpoints(controlFile,caseFile, d)
                        optT[i][folderList.index(f)][dilutionList.index(d)]=optt_s
                        cov[i][folderList.index(f)][dilutionList.index(d)]=cov_s
    
            cov=np.mean(np.mean(cov,0),1)
            with h5py.File(outputFile, 'w') as f:
                f.create_dataset('cov', data=cov)
                f.create_dataset('optT', data=optT)
                f.close()
                
        optT_avg=np.mean(optT,0)
        optT_std=np.std(optT,0)
    
        
        # plot the one dimensional figure 
        fig=plt.figure()
        for k in xrange(len(folderList)):
            ax=fig.add_subplot(2,2,k+1)
            ax.semilogx(dilutionList,optT_avg[k],'b*-',label=str(cov[k]))
            for i in xrange(len(dilutionList)):
                ax.semilogx([dilutionList[i],dilutionList[i]],
                [optT_avg[k][i]-optT_std[k][i],
                optT_avg[k][i]+optT_std[k][i]],
                color='r')
                
            if m==0 and k!=0:
                ax.legend(loc=3,title='read depth')
            else:
                ax.legend(loc=2,title='read depth')
                
            ax.set_xlabel('Dilution Rate')
            ax.set_ylabel('Optimal Threshold')
        title=figtitleopt[m]
       
        #plt.savefig(title)
        
        
        # plot the two dimensional figure
        # in log scale
        cov=[math.log(x)/math.log(10) for x in cov]
        dilutionList=[math.log(x)/math.log(10) for x in dilutionList]
        
        # spline

        f = interpolate.interp2d(dilutionList, cov, optT_avg, kind='cubic')
        steps=30
        covstep=float(cov[len(cov)-1]-cov[0])/steps      
        ncov = np.arange(cov[0], cov[len(cov)-1], covstep)
        
        dstep=float(dilutionList[len(dilutionList)-1]-dilutionList[0])/steps
        ndilutionList = np.arange(dilutionList[0], dilutionList[len(dilutionList)-1], dstep)
        noptT_avg = f(ndilutionList, ncov)
        
        # reshape the arrays
        rdilutionList,rcov =np.meshgrid(ndilutionList,ncov)      
        fig=plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_surface(rdilutionList, rcov, noptT_avg, rstride=1, cstride=1)
       
        ax.set_xlabel('log(dilution rate)')
        ax.set_ylabel('log(coverage)')
        ax.set_zlabel('optimal threshold')
              
        plt.savefig('%s_3D' %title)
        
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

