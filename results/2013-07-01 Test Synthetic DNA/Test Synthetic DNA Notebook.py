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

print os.getcwd()
import scipy.stats as ss

dilution=.3

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
    eeControl=f['ee'][...]

with h5py.File(caseFile, 'r') as f:
    phiCase = {}
    muCase_s = f['mu'][...]
    thetaCase_s = f['theta'][...]
    MCase = f['phi/M'][...]
    phiCase['mu0'] = f['phi/mu0'][()]
    phiCase['M0'] = f['phi/M0'][()]
    eeCase=f['ee'][...]
    

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


plt.figure()
plt.plot(SNRControl, SNRCase, marker='o', linestyle='None')
plt.plot(SNRControl[tpLoc-1], SNRCase[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Reference SNR')
plt.ylabel('Sample SNR')
plt.title('SNR Ratio for Variant Call')
plt.savefig("SNR Scatterplot.pdf")

## set the SNR shreshold as 1
threshold=1
indices = [i for i in range(J) if SNR[i]>1]
indices = np.array(indices)
tpLoc = np.arange(85,346,20)

# save data in a excel sheet
df = pd.DataFrame({'Postion':indices,'ee1': eeCase[1,indices,0],'ee2': eeCase[1,indices,1],'ee3': eeCase[1,indices,2]})
file_name='Dilution='+str(dilution)+' results.xlsx'
df.to_excel(file_name,sheet_name='sheet1', index=False)

iMutation=[]
for i in range(len(tpLoc)):
    idx=[idx for idx in range(len(indices)) if indices[idx]==(tpLoc[i]-1)]
    iMutation.append(idx)
iMutation=np.reshape(iMutation,len(iMutation))

lamda=2.0/3
m=len(indices)
p=np.zeros((m,N))
for i in range(m):
    for j in xrange(N):
        p[i,j]=rvd27.chi2test(eeCase[j,indices[i],:],lamda=lamda)

p=p+1e-30*np.ones_like(p)

p_logmax=np.log10(np.amax(p,axis=1))
p_logmin=np.log10(np.amin(p,axis=1))

plt.figure()
f,axarr = plt.subplots(5,sharex=True)
axarr[0].plot(roi, SNR, marker='o')
axarr[0].plot(tpLoc, SNR[tpLoc-1], color='r', marker='o', linestyle='None')
axarr[0].axhline(y=1, xmin=0, xmax=400,color='g',linestyle='--')
axarr[0].set_ylabel('SNR')
title='Dilution='+str(dilution)+'SNR Ratio and chi2 test for Variant Call'
axarr[0].set_title(title)

axarr[1].plot(indices, p_logmax, color='b', marker='o', linestyle='None')
axarr[1].plot(indices[iMutation], p_logmax[iMutation],color='r',marker='o',linestyle='None')
axarr[1].set_ylabel('log(pmax)')

axarr[2].plot(indices, p_logmin, color='b', marker='o', linestyle='None')
axarr[2].plot(indices[iMutation], p_logmin[iMutation],color='r',marker='o',linestyle='None')
axarr[2].set_ylabel('log(pmin)')

axarr[3].plot(indices, np.log10(p), color='b', marker='o', linestyle='None')
axarr[3].plot(indices[iMutation], np.log10(p[iMutation]),color='r',marker='o',linestyle='None')
axarr[3].set_ylabel('log(p)')

axarr[4].plot(indices, p, color='b', marker='o', linestyle='None')
axarr[4].plot(indices[iMutation], p[iMutation],color='r',marker='o',linestyle='None')
axarr[4].set_xlabel('Location')
axarr[4].set_ylabel('origin(p)')

file_name='Dilution='+str(dilution)+' SNR.pdf'
plt.savefig(file_name)
