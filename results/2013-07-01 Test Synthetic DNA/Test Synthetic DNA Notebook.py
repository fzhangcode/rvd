import sys
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import pickle
import multiprocessing as mp
import h5py
import logging
from scipy.special import gammainc
import pdb
print os.getcwd()
import scipy.stats as ss

dilution=0.3

if dilution==0.1:
    caseFile='case0_1.hdf5'
elif dilution==0.3:
    caseFile='case0_3.hdf5'
elif dilution==1:
    caseFile='case1_0.hdf5'
elif dilution==10:
    caseFile='case10_0.hdf5'
else:
    caseFile='case100_0.hdf5'

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

#os.chdir("/Users/pjflaherty/Research/rvd2/results/2013-07-01 Test Synthetic DNA")
#os.chdir("S:/yhe2/Research/rvd2/results/2013-07-01 Test Synthetic DNA")

tocfilename = "synthetic_toc.txt"
toc = pd.read_table(tocfilename)
##print toc

with h5py.File('control.hdf5', 'r') as f:
    phiControl = {}
    muControl_s = f['mu'][...]
    thetaControl_s = f['theta'][...]
    MControl = f['phi/M'][...]
    phiControl['mu0'] = f['phi/mu0'][()]
    phiControl['M0'] = f['phi/M0'][()]

with h5py.File(caseFile, 'r') as f:
    phiCase = {}
    muCase_s = f['mu'][...]
    thetaCase_s = f['theta'][...]
    MCase = f['phi/M'][...]
    phiCase['mu0'] = f['phi/mu0'][()]
    phiCase['M0'] = f['phi/M0'][()]

muControl = np.array(np.median(muControl_s, 1))
muCase = np.median(muCase_s,1)
mu0 = (phiCase['mu0'] + phiControl['mu0'])/2

(N,J,nsample) = np.shape(thetaCase_s)
roi = np.arange(J)
tpLoc = np.arange(85,346,20)


plt.figure()
plt.plot(roi, muCase/muControl, marker='o')
plt.plot(tpLoc, muCase[tpLoc-1]/muControl[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Location')
plt.ylabel('mu Case / mu Control')
plt.title('Case/Control Ratio for Variant Call')
plt.savefig("case-control.pdf")

iqrControl = np.percentile(muControl_s, (2.5, 97.5), axis=1)
iqrCase = np.percentile(muCase_s, (2.5, 97.5), axis=1)

SNRControl = muControl/(iqrControl[1] - iqrControl[0] + 1/np.sqrt(phiControl['M0']))
SNRCase = muCase/(iqrCase[1] - iqrCase[0] + 1/np.sqrt(phiCase['M0']))
SNR = SNRCase/SNRControl

plt.figure()
plt.plot(roi, SNR, marker='o')
plt.plot(tpLoc, SNR[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Location')
plt.ylabel('SNR')
plt.title('SNR Ratio for Variant Call')
plt.savefig("SNR.pdf")

plt.figure()
plt.plot(SNRControl, SNRCase, marker='o', linestyle='None')
plt.plot(SNRControl[tpLoc-1], SNRCase[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Reference SNR')
plt.ylabel('Sample SNR')
plt.title('SNR Ratio for Variant Call')
plt.savefig("SNR Scatterplot.pdf")

def load_sequencing(dcFileNameList):
    r=[]; n=[];e=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    for dcFileName in dcFileNameList:
        with open(dcFileName, 'r') as dcFile:
            header = dcFile.readline().strip()
            dc = dcFile.readlines()
            dc = [x.strip().split("\t") for x in dc]
            loc = map(int, [x[2] for x in dc])
            refb = [x[4] for x in dc]
            c = [map(int, x[5:9]) for x in dc]
            c = np.array(c)
            (J, K) = np.shape(c)
            n1 = np.sum(c, 1)
            r1 = np.zeros(J)
            Inx=np.zeros(J)
            for j in xrange(0, J):
                r1[j] = n1[j] - c[j, acgt[refb[j]]]
                Inx[j] = 4*j+acgt[refb[j]]
            c=np.delete(c,Inx,None)
            c=np.reshape(c,(-1,3))    
        r.append(r1)
        n.append(n1)
        e.append(c)
    r = np.array(r)
    n = np.array(n)
    e = np.array(e)
    return (r,n,e)

logging.debug("Read in sequencing depth of dilution: %0.1f" % dilution)
caseFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
(r, n, e) = load_sequencing(caseFileList)

print np.shape(e)
print np.shape(r)
print np.shape(n)



def chisquaretest(X, lamda=2.0/3, pvector=np.array([1.0/3]*3)):
    nsum=np.sum(X)
    E=nsum*pvector
    X=np.array(X)

    if lamda==0 or lamda==-1:
        C=2.0/np.sum(X*np.log(X*1.0/E))
    else:
        C=2.0/(lamda*(lamda+1))*np.sum(X*((X*1.0/E)**lamda-1))
    df=len(pvector)-1
    #p=scipy.special.gammainc(C,df)
    # p=1-gammainc(df/2,C/2)
    p = 1 - ss.chi2.cdf(C, df) 
    return(p)


## set the SNR shreshold as 1
threshold=1
indices = [i for i in range(J) if SNR[i]>1]


m=len(indices)
p=np.zeros((m,N))
for i in range(m):
    for j in range(N):
        p[i,j]=chisquaretest(e[j,indices[i],:])

print e[0,indices,:]

print indices
print p       
pdb.set_trace()
pmedian=np.median(p,1)


##plt.figure()
##plt.plot(indices, pmedian, marker='o')
##plt.plot(tpLoc, SNR[tpLoc-1], color='r', marker='o', linestyle='None')
##plt.xlabel('Location')
##plt.ylabel('SNR')
##plt.title('SNR Ratio for Variant Call')
##plt.savefig("SNR.pdf")

