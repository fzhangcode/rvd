import sys
import os
import numpy as np
import h5py
import logging
import pdb
import pandas as pd
import random

os.chdir('S:/yhe2/Research/rvd2/results/2013-11-07_HCC1187_PAXIP1_genome_fa')
# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../../rvd2/src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    random.seed(199096)
    N=1000 # Z sampling size  
            
    # Load the control model samples
    controlFile="Control.hdf5"
    (controlPhi, controlTheta, controlMu, controlloc, controlR, controlN) = rvd27.load_model(controlFile)
    
    caseFile="Case.hdf5"
    (casePhi, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd27.load_model(caseFile)
    
    # Extract the common locations in case and control
    caseLocIdx = [i for i in xrange(len(caseLoc)) if caseLoc[i] in controlloc]
    controllocIdx = [i for i in xrange(len(controlloc)) if controlloc[i] in caseLoc]

    caseMu = caseMu[caseLocIdx,:]
    controlMu = controlMu[controllocIdx,:]
    caseR = caseR[:,caseLocIdx,:]
    controlR = controlR[:,controllocIdx,:]
    loc = caseLoc[caseLocIdx]
    J = len(loc)

    # caculation the probability that mu is in range (0, 0.1),(0.05,0.9),(0.9,1) respecitively.
    somatic_interval=[(0, 0.1),(0.05,0.9),(0.9,1)]
    muControlP = Mu_distr(controlMu-controlPhi['mu0'], somatic_interval)
    muCaseP = Mu_distr(caseMu-casePhi['mu0'], somatic_interval)

    # For each position find the range with highest probability
    ControlMaxpIdx = np.argmax(muControlP, axis=1)
    CaseMaxpIdx = np.argmax(muCaseP,axis=1)

    pdb.set_trace()
    # Germline mutation classify; only looking into controlmu
    Germline_type={
        0:'Reference',
        1:'Heterozygous',
        2:'Homozygous',
        }
    Germline_status=[]
    for i in xrange(J):
        Germline_status.append(Germline_type[ControlMaxpIdx[i]])
    Germline_status=np.array(Germline_status)
    

    # Somatic mutation classify
    Somatic_type={
        (0,0):'Reference',
        (0,1):'Heterozygous',
        (0,2):'Homozygous',
        (1,0):'LOH',
        (1,1):'Unknown',
        (1,2):'LOH',
        (2,0):'Homozygous',
        (2,1):'Heterozygous',
        (2,2):'Unknown',
        }
    Somatic_status=[]
    for i in xrange(J):
        Somatic_status.append(Somatic_type[(ControlMaxpIdx[i],CaseMaxpIdx[i])])
    Somatic_status=np.array(Somatic_status)

    # non-reference index
    nonref_idx=[]
    for i in xrange(J):
        #if postP[i]>=0.95 and chi2P[i]<0.05/J:
        #if postP[i]>=0.95:
        if Germline_status[i]!='Reference' or Somatic_status[i]!='Reference':
            nonref_idx.append(i)
    nonref_idx=np.array(nonref_idx)


    hdf5file='postTest.hdf5'
    try:
        with h5py.File('postTest.hdf5','r') as f:
            postP=f['postP'][...]
            T=f['T'][()]
            chi2P=f['chi2pvalue'][...]
            f.close()
    except IOError:     
        T = 0.00001 # detection threshold
        print("Detection Threshold = %f" % T)
        
        # Sample from the posterior Z = muCase - muControl
        (Z, caseMuS, controlMuS) = rvd27.sample_post_diff(caseMu, controlMu, N)
        # Bayesian Hypothesis Testing
        postP = rvd27.bayes_test(Z, [(T, np.inf)]) # Posterior Prob that muCase is greater than muControl by T

        # Chi2 Testing
        nRep = caseR.shape[0]
        chi2P = np.zeros((J,nRep))
        for j in xrange(J):
                chi2P[j,:] = np.array([rvd27.chi2test( caseR[i,j,:] ) for i in xrange(nRep)] )
        # Save the test results
        f = h5py.File(hdf5file)
        f.create_dataset('loc', data=loc)
        f.create_dataset('postP', data=postP)
        f.create_dataset('T', data=T)
        f.create_dataset('chi2pvalue',data=chi2P)
        f.close()


    pdb.set_trace()
    print('Done')
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

    idx=nonref_idx
    with open('result.txt','w') as f:
        print >> f, 'count\tloc\tpostP\tGermline_Status\tSomatic_Status\trefb\tA0\tC0\tG0\tT0\tA1\tC1\tG1\tT1'
        if len(idx):
            for f0, f1, f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13 in zip(idx,loc[idx],postP[idx],
                                                                   Germline_status[idx],
                                                                   Somatic_status[idx],
                                                                   refb[idx],
                                                                   A0[idx],C0[idx],
                                                                   G0[idx],T0[idx],
                                                                   A1[idx],C1[idx],
                                                                   G1[idx],T1[idx],):
                print >> f, f0,'\t', f1,'\t',f2,'\t',f3,'\t',f4,'\t',f5,'\t',f6,\
                      '\t',f7,'\t',f8,'\t',f9,'\t',f10,'\t',f11,'\t',f12,'\t',f13
        f.close()
        

def Mu_distr(Z, roi):
        """ Return posterior probabilities in regions defined in list of tuples (roi)
            from samples in columns of Z. """

        (J,N)=np.shape(Z)
        
        nTest = len(roi) # get the number of regions to compute probabilities 
        
        p = np.zeros((J,nTest))
        for i in xrange(nTest):
            for j in xrange(J):
                    p[j,i] = np.float( np.sum( np.logical_and( (Z[j,:] >= roi[i][0]),\
                                                               (Z[j,:] <= roi[i][1]) ) ) ) / N

        return p

if __name__ == "__main__":
	main()
