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
import scipy
import pdb
print os.getcwd()
import scipy.stats as ss

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

with h5py.File('control.hdf5', 'r') as f:
    muControl_s = f['mu'][...]
(J,N)=np.shape(muControl_s)


dilution_opt=(0.1,0.3,1.0,10.0,100.0)
t=[]
var=False # whether the variance of muCase_s and muControl_s is assumed equal or not
for d in xrange(5):
    dilution=dilution_opt[d]
    
    if dilution==0.1:
        caseFile='case0_1.hdf5'
    elif dilution==0.3:
        caseFile='case0_3.hdf5'
    elif dilution==1.0:
        caseFile='case1_0.hdf5'
    elif dilution==10.0:
        caseFile='case10_0.hdf5'
    else:
        caseFile='case100_0.hdf5'

    with h5py.File(caseFile, 'r') as f:
        muCase_s = f['mu'][...]
    
    ## plot histogram
    position1=49
    fig=plt.figure()
    ax1=fig.add_subplot(2,1,1)
    ax1.hist(muCase_s[position1,:],20)
    ax1.hist(muControl_s[position1,:],20)
    title='Histogram of mu when '+'dilution='+str(dilution)+' position='+str(position1+1)
    ax1.set_title(title)
    ax1.set_ylabel('t')
    ax1.legend( ['Case','Control'])

    position2=244
    ax2=fig.add_subplot(2,1,2)
    ax2.hist(muCase_s[position2,:],20)
    ax2.hist(muControl_s[position2,:],20)
    title='position='+str(position2+1)
    ax2.set_title(title)
    ax2.set_ylabel('t')
    ax2.set_xlabel('mu')
    ax2.legend(['Case','Control'])
    plt.savefig('Histogram of mu_s when Dilution='+str(dilution)+'.jpg')
    

    ##    scipy.stats.ttest_ind(a, b, axis=0, equal_var=True)[source] returns
    ##    a test decision for the null hypothesis that the data in vectors xand
    ##    y comes from independent random samples normal distributions with
    ##    equal means and equal but variances,using the two-sample t-test.
    stat=np.array(scipy.stats.ttest_ind(muCase_s, muControl_s, axis=1, equal_var=var))
    t.append(stat[0,:])
    
t=np.array(t)    
pos=np.arange(J)
mpos=np.arange(85,346,20)

fig2=plt.figure(figsize=(9,16))
for i in xrange(5):
    ax=fig2.add_subplot(3,2,i+1)
    ax.plot(pos,t[i,:],marker='o')
    ax.plot(mpos-1,t[i,mpos-1],color='r',marker='o',linestyle='None')
    ax.set_xlabel('Location')
    ax.set_ylabel('t')
    ax.set_title('Dilution='+str(dilution_opt[i]))
    
if var:
    title='tstat '+'with assumed equal variace'+'.jpg'
else:
    title='tstat '+'with assumed unequal variace'+'.jpg'
    
plt.savefig(title)
