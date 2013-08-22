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
    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)

    Nreads=[]
    label=[]       
    for d in dilutionList:
        
        h5FileName = 'ROC_dilution%d.hdf5' % d
        logging.debug("Processing dilution: %0.1f%%" % d)        

        controlFile="Control.hdf5"
        caseFile="Case%s.hdf5" % str(d).replace(".","_")

        [fpr, tpr, reads] = ROCpoints(controlFile,caseFile)
        Nreads.extend(reads)

        # ROC 
        plt.plot(fpr,tpr)
        label.extend(('%0.1f%%' % d,))
        
        
##    plt.title('Median of read depth=%d' % np.median(np.array(Nreads),axis=None))
    plt.legend(label,loc=4)
    plt.plot([0,1],[0,1],color='k',linestyle='dashed')

    plt.title('Read depth median = %d' % np.median(np.array(Nreads),axis=None))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.savefig('ROC_plot.png')

 
##        """ plot the histograms"""
##        position=[8,65,106,387,205,225]
##        title = "histZ_dilution=%s" % str(d).replace(".","_")
##        histPlotZ(position,Z,title,figform='.png')
##                
##                
##        title='histXY_dilution=%s' % str(d).replace(".","_")
##        histPlotXY(position, caseMuS, controlMuS, title, figform='.png')
#                ROCsubplot(figroc,points,d+1,'dilution='+str(dilution))
#            figroc.savefig(titleroc+figform)        

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
    pos = np.arange(85,346,20)
    pos = [loc[i] for i in xrange(J) if loc[i] in pos]

    
    # Sample from the posterior Z = muCase - muControl        
    (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu, controlMu, N)
    
    # Compute cumulative posterior probability for regions (Threshold,np.inf) 
    pList = [rvd27.bayes_test(Z, [(T, np.inf)]) for T in np.linspace(np.min(np.min(Z)), np.max(np.max(Z)), num=300)]

    # mutation classification
    clsList = np.array((np.array(pList)>P).astype(int))
    clsList = clsList.reshape((clsList.shape[0],clsList.shape[1]))# category list

    # false postive rate
    fpr = np.array([float(sum(clsList[i])-sum(clsList[i,np.array(pos)-1]))/(clsList.shape[1]-len(pos)) \
           for i in xrange(clsList.shape[0])])
    
    # true positive rate
    tpr = np.array([float(sum(clsList[i,np.array(pos)-1]))/len(pos) for i in xrange(clsList.shape[0])])

    reads=[caseN,controlN]
    return fpr,tpr,reads
        
    
def histPlotZ(position,Z,title=None,bins=40,subplotsize=[3,2],figform='.pdf'):
    fighandle=plt.figure(figsize=(16,9))
    plt.suptitle(title)
    for i in xrange(len(position)):        
        p=position[i]
        ax=fighandle.add_subplot(subplotsize[0],subplotsize[1],i+1)
        ax.hist(Z[p-1,:],bins)
        subtitle='position='+str(p)
        ax.set_title(subtitle)
    
    title=title+figform
    plt.savefig(title)

        
def histPlotXY(position, 
                muCase_s, 
                muControl_s, 
                title=None, 
                bins=40,
                subplotsize=[3,2],
                figform='.pdf'):
                
    fighandle=plt.figure(figsize=(16,9))
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


    '''ROC plot'''    
def ROCsubplot(figroc,points,NOsubplot,subtitle):
    ax=figroc.add_subplot(3,2,NOsubplot)
    ax.plot(points[:,0],points[:,1],color='b',marker='o')

    ax.set_title(subtitle)
    ax.grid(True)
    ax.set_ylabel('True Positive Rate')
    ax.set_xlabel('False Positive Rate')

    '''text the threshold when the ROC point is closest to the (0,1)point'''    
    dpoints=[points[:,0],points[:,1]-1]
    dpoints=np.power(dpoints,2)
    distance=np.sum(dpoints,0)
    ind=distance.argmin()
    ax.plot(points[ind,0],points[ind,1],color='r',marker='o')
    ax.text(0.2,0.6,'TH='+str('%0.4f'% points[ind,2]))
    ax.text(0.2,0.5,'ACC='+str('%0.4f'% points[ind,3]))
    ax.text(0.2,0.4,'FPR='+str('%0.4f'% points[ind,0]))
    ax.text(0.2,0.3,'TPR='+str('%0.4f'% points[ind,1]))

def save_rocPoints(h5FileName,muCase_s=None,muControl_s=None,Z=None,points=None,steps=None,alpha=None):

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
    if steps is not None:
        h5file.create_dataset('steps', data=steps)
    if alpha is not None:
        h5file.create_dataset('alpha', data=alpha)
    
    h5file.close()



if __name__ == '__main__':
    main()

