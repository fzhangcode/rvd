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
    dilutionList = (0.1,0.3,1.0,10.0)

    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg10',\
                  '2013-08-14_Compute_ROC_Synthetic_avg100',\
                  '2013-08-14_Compute_ROC_Synthetic_avg1000',\
                  '2013-08-14_Compute_ROC_Synthetic_avg10000',\
                  '2013-08-14_Compute_ROC_Synthetic_avg100000')
##    dilutionList = (1.0,)
##
##    folderList = ('2013-08-14_Compute_ROC_Synthetic_avg100000',)

    
    N=1000 # Z sampling size  
    fig=plt.figure(figsize=(10, 9))
    chi2=True
    plt.suptitle('Comparing ROC plot for different read depth' )
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        ax = fig.add_subplot(2,2,dilutionList.index(d)+1)
        label=[]
        for f in folderList:            
            controlFile = "../%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "../%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            
            [fpr,tpr, cov, optimalT, CalledCaseMu, CalledControlMu, CalledLoc, predList, trueList] = ROCpoints(controlFile,caseFile,P=0.95,chi2=chi2)
            title='MuBarPlot_Dilution%(dilution)0.1f_Depth%(depth)d.png'%{'dilution':d,'depth':cov}
            
            #MuBarPlot
            MuBarPlot(CalledControlMu,CalledCaseMu,CalledLoc,optimalT,title)
            
            #Confusion Matrix
            title='ConfusionMatrix_Dilution%(dilution)0.1f_Depth%(depth)d.png'%{'dilution':d,'depth':cov}
            Conf_Matrix(predList,trueList,title)

            # ROC curve
            ax.plot(fpr,tpr, label='%d' % cov)

        l = ax.legend(loc=4,prop={'size':9},title='Read depth')
        l.get_title().set_fontsize(9)
        ax.plot([0,1],[0,1],color='k',linestyle='dashed')
        ax.set_title('%0.1f%% Mutant Mixture' % d)
        ax.set_xlim((-0.01,1.01))
        ax.set_ylim((-0.01,1.01))
                 
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        title = 'Dilution%0.1f%%' % d
        ax.set_title(title, fontsize=10)
        
    plt.sca(ax)
    if chi2:
        title='ROC_with_chi2'
    else:
        title='ROC_without_chi2'
    plt.savefig(title)
    
def ROCpoints(controlFile,caseFile,N=1000,P=0.95,chi2=False):
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
    optimalT = T[Tidx]
    CalledCaseMu = caseMu[clsList[Tidx]==1]
    CalledControlMu = controlMu[clsList[Tidx]==1]
    CalledLoc = loc[clsList[Tidx]==1]

    # return information for confusion matrix at optimal threshold optimalT
    
    
    predList = clsList[Tidx]# prediction classification
    trueList = np.zeros_like(predList)# true classification
    trueList[posidx] =np.ones_like(posidx)
    return fpr,tpr, cov, optimalT, CalledCaseMu, CalledControlMu, CalledLoc, predList, trueList

def MuBarPlot(CallControlMu,CallCaseMu,Loc,optimalT,title):
    if len(Loc) is not 0:
        ind = np.arange(CallCaseMu.shape[0])
        width=0.35

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, np.mean(CallControlMu,1), width, color='0.8', yerr=np.std(CallControlMu,1))
        rects2 = ax.bar(ind+width, np.mean(CallCaseMu,1), width, color='0.5', yerr=np.std(CallCaseMu,1))

        # add some
        ax.set_ylabel('Minor Allele Frequency')
        ax.set_xlabel('Called Locations')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( [str(x)[8:] for x in Loc] )
        ax.set_ylim
        ax.legend( (rects1[0], rects2[0]), ('Control', 'Case') )
        fig.suptitle('Optimal Threshold = %0.5f'%optimalT )
        plt.savefig(title.replace('.','_',1))
    else:
        print 'No variant is called'
        
def Conf_Matrix(predList,trueList,title):
    L = len(predList)
    TN = len([i for i in xrange(L) if predList[i]==0 and trueList[i]==0])
    FP = len([i for i in xrange(L) if predList[i]==1 and trueList[i]==0])
    FN = len([i for i in xrange(L) if predList[i]==0 and trueList[i]==1])
    TP = len([i for i in xrange(L) if predList[i]==1 and trueList[i]==1])
##    pdb.set_trace()
    conf_arr=[[TN,FP],[FN,TP]]
    norm_conf = []
    for i in conf_arr:
            a = 0
            tmp_arr = []
            a = sum(i,0)
            for j in i:
                    tmp_arr.append(float(j)/float(a))
            norm_conf.append(tmp_arr)
    label_conf=[['TNR','FPR'],['FNR','TPR']]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    res = ax.imshow(np.array(norm_conf), cmap=plt.cm.jet, interpolation='nearest')
    for i,j in ((x,y) for x in xrange(len(conf_arr))
                for y in xrange(len(conf_arr[0]))):
        ax.annotate(str(conf_arr[i][j]),xy=(i,j-0.1),color='white')
        label = '%(l)s = %(r)0.2f%%' %{'l':label_conf[i][j],'r':norm_conf[i][j]*100}
        ax.annotate(label,xy=(i-0.15,j+0.1),color='white')

    cb = fig.colorbar(res)
    plt.title('Confusion matrix')
    plt.xticks(np.arange(0, 2, 1))
    plt.yticks(np.arange(0, 2, 1))
    plt.xlabel('True outcome')
    plt.ylabel('Prediction outcome')
    plt.savefig(title.replace('.','_',1))
if __name__ == '__main__':
    main()
