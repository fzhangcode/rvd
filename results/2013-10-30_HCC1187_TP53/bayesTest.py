import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import pandas as pd

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../../rvd2/src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    N=1000 # Z sampling size  
            
    # Load the control model samples
    controlFile="Control.hdf5"
    (controlPhi, controlTheta, controlMu, controlloc, controlR, controlN) = rvd27.load_model(controlFile)
    
    caseFile="Case.hdf5"
    (casePhi, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd27.load_model(caseFile)
    
    T = 0 # detection threshold
    print("Detection Threshold = %f" % T)

    # Extract the common locations in case and control
    caseLocIdx = [i for i in xrange(len(caseLoc)) if caseLoc[i] in controlloc]
    controllocIdx = [i for i in xrange(len(controlloc)) if controlloc[i] in caseLoc]

    caseMu = caseMu[caseLocIdx,:]
    controlMu = controlMu[controllocIdx,:]
    caseR = caseR[:,caseLocIdx,:]
    controlR = controlR[:,controllocIdx,:]
    loc = caseLoc[caseLocIdx]
    J = len(loc)

    # Sample from the posterior Z = muCase - muControl        
    (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu, controlMu, N)

    postP = rvd27.bayes_test(Z, [(T, np.inf)]) # Posterior Prob that muCase is greater than muControl by T
	
    nRep = caseR.shape[0]
    chi2P = np.zeros((J,nRep))
    for j in xrange(J):
	    chi2P[j,:] = np.array([rvd27.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )

    hdf5file='postTest.hdf5'
    try:
        os.remove(hdf5file)
    except OSError:
        pass
    
    # Save the test results
    f = h5py.File(hdf5file)
    f.create_dataset('loc', data=loc)
    f.create_dataset('postP', data=postP)
    f.create_dataset('T', data=T)
    f.create_dataset('chi2pvalue',data=chi2P)
    f.close()

    idx = [i for i, j in enumerate(postP) if j>=0.95]

    # Read in called position depth chart
    contdc=pd.read_table('depth_chart/HCC1187BL_S1.dc')
    casedc=pd.read_table('depth_chart/HCC1187C_S1.dc')

    NonN=[i for i, j in enumerate(contdc.refb) if j is not 'N']
    
    refb=np.array(contdc.refb[NonN])
    A0=np.array(contdc.A_All[NonN])
    C0=np.array(contdc.C_All[NonN])
    G0=np.array(contdc.G_All[NonN])
    T0=np.array(contdc.T_All[NonN])

    A1=np.array(casedc.A_All[NonN])
    C1=np.array(casedc.C_All[NonN])
    G1=np.array(casedc.G_All[NonN])
    T1=np.array(casedc.T_All[NonN])

    with open ('controlR.txt','w') as f:
        for i in xrange(np.shape(controlR)[1]):
            print >> f, i+1,'\t',controlR[0,i,0],'\t',controlR[0,i,1],'\t',controlR[0,i,2]
        f.close()

    
    with open ('caseR.txt','w') as f:
        for i in xrange(np.shape(caseR)[1]):
            print >> f, i+1,'\t', caseR[0,i,0],'\t',caseR[0,i,1],'\t',caseR[0,i,2]
        f.close()   
    with open('result.txt','w') as f:
        print >> f, 'count\tloc\tpostP\trefb\tA0\tC0\tG0\tT0\tA1\tC1\tG1\tT1'
        for f0, f1, f2,f3,f4,f5,f6,f7,f8,f9,f10,f11 in zip(idx,loc[idx], postP[idx],refb[idx],
                                                       A0[idx],C0[idx],
                                                       G0[idx],T0[idx],
                                                       A1[idx],C1[idx],
                                                       G1[idx],T1[idx],):
            print >> f, f0,'\t', f1,'\t',f2,'\t',f3,'\t',f4,'\t',f5,'\t',f6,'\t',f7,'\t',f8,'\t',f9,'\t',f10,'\t',f11
        f.close()
        
if __name__ == "__main__":
	main()
