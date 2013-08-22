import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    N=1000 # Z sampling size  
            
    # Load the control model samples
    controlFile="Control.hdf5"
    (controlPhi, controlTheta, controlMu, controlLoc, controlR, controlN) = rvd27.load_model(controlFile)
    
    T = 0.0005 # detection threshold
    print("Detection Threshold = %f" % T)
    for d in dilutionList:
        
        print("Processing dilution: %0.1f%%" % d)
        
        # Load the case model samples
        caseFile="Case%s.hdf5" % str(d).replace(".","_")
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

        # Sample from the posterior Z = muCase - muControl        
        (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu, controlMu, N)

        postP = rvd27.bayes_test(Z, [(T, np.inf)]) # Posterior Prob that muCase is greater than muControl by T
	
        nRep = caseR.shape[0]
	chi2P = np.zeros((J,nRep))
	for j in xrange(J):
	    chi2P[j,:] = np.array([rvd27.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
	# Save the test results
	f = h5py.File('postTest%s.hdf5' % str(d).replace('.','_'),'w')
	f.create_dataset('loc', data=loc)
	f.create_dataset('postP', data=postP)
	f.create_dataset('T', data=T)
	f.create_dataset('chi2pvalue',data=chi2P)
	f.close()


        # Print the results for a select group of positions
        (nRep, J, M) = caseR.shape
	testPos = (8, 65, 106, 387, 205, 225)
        for pos in testPos:
            j = pos-1
            chi2P = [rvd27.chi2test( caseR[i,j,:] ) for i in xrange(nRep)]
            print("Position %d: Mutation Probability=%0.2f" % (caseLoc[j], postP[j]) )
	    print "Mismatch Base Uniformity Test p-value:",
	    for i in xrange(nRep):
		    print " %0.2f" % chi2P[i],
            print

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

