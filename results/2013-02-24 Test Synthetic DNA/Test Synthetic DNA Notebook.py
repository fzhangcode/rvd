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

# Insert the src/python/rvd26 directory at front of the path
rvddir = os.path.join('../../src/python/rvd26')
sys.path.insert(0, rvddir)
import rvd26

tocfilename = "synthetic_toc.txt"
toc = pd.read_table(tocfilename)
print toc

with h5py.File('control.hdf5', 'r') as f:
    phiControl = {}
    muControl_s = f['mu_s'][...]
    thetaControl_s = f['theta_s'][...]
    MControl_s = f['M_s'][...]
    phiControl['a'] = f['phi/a'][()]
    phiControl['b'] = f['phi/b'][()]
    phiControl['mu0'] = f['phi/mu0'][()]
    phiControl['M0'] = f['phi/M0'][()]

muControl = np.array(np.median(muControl_s, 1))
muCase = np.median(muCase_s,1)
mu0 = (phiCase['mu0'] + phiControl['mu0'])/2

(N,J,nsample) = np.shape(thetaCase_s)
roi = np.arange(J)
tpLoc = np.arange(85,346,20)


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

plt.plot(roi, SNR, marker='o')
plt.plot(tpLoc, SNR[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Location')
plt.ylabel('SNR')
plt.title('SNR Ratio for Variant Call')
plt.savefig("SNR.pdf")

plt.plot(SNRControl, SNRCase, marker='o', linestyle='None')
plt.plot(SNRControl[tpLoc-1], SNRCase[tpLoc-1], color='r', marker='o', linestyle='None')
plt.xlabel('Reference SNR')
plt.ylabel('Sample SNR')
plt.title('SNR Ratio for Variant Call')
plt.savefig("SNR Scatterplot.pdf")