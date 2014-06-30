# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import pdb

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3, venn2
import vcf
import re
import h5py

def main():
    # read in the vcf files
    filename = 'T0diploid_S2_0_01.vcf'
    t = pd.read_table(filename, header = 1, skiprows = 9)
    T0POS = t.POS
    T0S1AF = map(AFfun, t.INFO)
    # S1AF = [float(s) for s in infonum]

    ###################### plot the MAF levels in control and case sample for each experiment.########################  
    ## T1S3 diploid case; T0S2 diploid control
    filename = 'T1diploid_S3.hdf5'
    T0S2AF, T1S3AF, T1T0POS= hdf5_import(filename)
    bar2(T0S2AF, T1S3AF, T1T0POS, r = 25, slegend = ['control T0S2', 'case T1'], stitle ='T1_T0S2.png')

    filename = 'T2diploid_S4.hdf5'
    T0S2AF_2, T2S4AF, T2T0POS= hdf5_import(filename)
    bar2(T0S2AF_2, T2S4AF, T2T0POS, r = 45, slegend = ['control T0S2', 'case T2'], stitle ='T2_T0S2.png')

    filename = 'T1diploid_S3_T2diploid_S4.hdf5'
    S3AF, S4AF, T2T1POS= hdf5_import(filename)
    bar2(S3AF, S4AF, T2T1POS, r = 25, slegend = ['control T1', 'case T2'], stitle ='T2_T1.png')

    filename = 'T0haploid_S1.hdf5'
    S2AF, S1AF, S2S1POS= hdf5_import(filename)
    bar2(S2AF, S1AF, S2S1POS, r = 25, slegend = ['control T0S2', 'case T0S1'], stitle ='T0S2_T0S1.png')

    ########################################### Plot the positions as Venn diagram #####################################  
    plt.figure(1)
    plt.subplot(212)
    venn3([set(T1T0POS), set(T2T0POS),set(T2T1POS)], ('T1-T0_S2','T2-T0_S2','T2-T1'))
    plt.subplot(211)
    venn3([set(T0POS), set(T1T0POS), set(T2T0POS)], ('T0_S2','T1-T0_S2','T2-T0_S2'))
    plt.savefig('vennfig.png')


    #################### Plot the common positions for (T0 diploid control, T1 and T2 case) as bar plot#################
    T0T1T3POS = list(set(T1T0POS) & set(T2T0POS))
    T0T1T3POS = sorted(T0T1T3POS)
    L = len(T0T1T3POS)
    T0AF = np.array([T0S2AF[list(T1T0POS).index(T0T1T3POS[i])] for i in xrange(L)])
    T1AF = np.array([T1S3AF[list(T1T0POS).index(T0T1T3POS[i])] for i in xrange(L)])
    T2AF = np.array([T2S4AF[list(T2T0POS).index(T0T1T3POS[i])] for i in xrange(L)])

    bar3(T0AF, T1AF, T2AF, T0T1T3POS, r = 25, slegend = ['T0S2', 'T1', 'T2'], stitle ='T0S2T1T3POSbar.png')

def bar3(MAF1, MAF2, MAF3, pos, r = 25, slegend = ['T0S2', 'T1', 'T2'], stitle ='MAF3.png'):
    alpha = 0.05
    cred = int(MAF1.shape[1]*alpha/2)
    #pdb.set_trace()
    # sort along Gibbs samples
    sortMAF1 = np.sort(MAF1,axis=1)
    yerr1 = np.array([np.mean(MAF1,1)-sortMAF1[:,cred], sortMAF1[:,MAF1.shape[1]-cred]-np.mean(MAF1,1)])

    sortMAF2 = np.sort(MAF2,axis=1)
    yerr2 = np.array([np.mean(MAF2,1)-sortMAF2[:,cred], sortMAF2[:,MAF2.shape[1]-cred]-np.mean(MAF2,1)])    

    sortMAF3 = np.sort(MAF3,axis=1)
    yerr3 = np.array([np.mean(MAF3,1)-sortMAF3[:,cred], sortMAF3[:,MAF3.shape[1]-cred]-np.mean(MAF3,1)])  

    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)  
    L = len(pos)
    ind = np.arange(L)
    width = 0.5
    rects1 = ax.bar(ind, np.mean(MAF1,1), width/2.0,
                    color='k', yerr=yerr1,ecolor='k')
    rects2 = ax.bar(ind+width/2.0, np.mean(MAF2,1), width/2.0,
                        color='r', yerr=yerr2,ecolor='k')
    rects2 = ax.bar(ind+width, np.mean(MAF3,1), width/2.0,
                        color='b', yerr=yerr3,ecolor='k')
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel('MAF (%)')
    ax.set_xlabel('Positions ChrXV:780000+X')
    # pdb.set_trace()
    xTickMarks = [str(x-780000) for x in pos]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=r, fontsize=10)
    ax.legend( (rects1, rects2, rects3), slegend,loc=2 ,prop={'size':10})
    plt.tight_layout()    
    plt.savefig(stitle)



def bar2( MAF1, MAF2, pos, r = 25, slegend = ['control', 'case'], stitle ='MAF.png'):
    ## compute the 95% credible interval
    alpha = 0.05
    cred = int(MAF1.shape[1]*alpha/2)
    #pdb.set_trace()
    # sort along Gibbs samples
    sortMAF1 = np.sort(MAF1,axis=1)
    yerr1 = np.array([np.mean(MAF1,1)-sortMAF1[:,cred], sortMAF1[:,MAF1.shape[1]-cred]-np.mean(MAF1,1)])

    sortMAF2 = np.sort(MAF2,axis=1)
    yerr2 = np.array([np.mean(MAF2,1)-sortMAF2[:,cred], sortMAF2[:,MAF2.shape[1]-cred]-np.mean(MAF2,1)])    

    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)  
    L = len(pos)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax.bar(ind, np.mean(MAF1,1), width,
                    color='k', yerr=yerr1,ecolor='k')
    rects2 = ax.bar(ind+width, np.mean(MAF2,1), width,
                        color='r', yerr=yerr2,ecolor='k')
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel('MAF (%)')
    ax.set_xlabel('Positions ChrXV:780000+X')
    # pdb.set_trace()
    xTickMarks = [str(x-780000) for x in pos]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=r, fontsize=10)
    ax.legend( (rects1, rects2), slegend,loc=2 ,prop={'size':10})
    plt.tight_layout()    
    plt.savefig(stitle)


def hdf5_import(filename):
    with h5py.File(filename,'r') as f:
        call = f['call'][...]
        caseMu = f['caseMu'][...]
        controlMu = f['controlMu'][...]
        loc = f['loc'][...]
    caseMu = caseMu[call]
    controlMu = controlMu[call]
    loc = loc[call]
    loc = [int(x.split(':')[1]) for x in loc]

    return controlMu, caseMu, loc

def table_import(filename):
    t = pd.read_table(filename, header = 1, skiprows = 10)
    AF = map(AFfun, t.INFO)
    AF = np.array(AF) 

    return t.POS, AF

def AFfun(infostr):
    infonum = re.findall(r'[+-]?\d+.\d+', infostr)
    AFnum = [float(s) for s in infonum]
    return AFnum

if __name__ == "__main__":
    main()