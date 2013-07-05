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

# Insert the src/python/directory at front of the path
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
    eeControl['ee']=f['ee'][...]

with h5py.File(caseFile, 'r') as f:
    phiCase = {}
    muCase_s = f['mu'][...]
    thetaCase_s = f['theta'][...]
    MCase = f['phi/M'][...]
    phiCase['mu0'] = f['phi/mu0'][()]
    phiCase['M0'] = f['phi/M0'][()]
    eeCase['ee']=f['ee'][...]
    

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


##logging.debug("Read in sequencing depth of dilution: %0.1f" % dilution)
##caseFileList = ["../../data/synthetic_dcs/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
##(r, n, e) = runall.load_depth(caseFileList)

print np.shape(eeCase)
print np.shape(r)
print np.shape(n)



## set the SNR shreshold as 1
threshold=1
indices = [i for i in range(J) if SNR[i]>1]


m=len(indices)
p=np.zeros((m,N))
for i in range(m):
    for j in range(N):
        p[i,j]=rvd27.chi2test(eeCase[j,indices[i],:])

print eeCase[0,indices,:]

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

