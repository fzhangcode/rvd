# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# import sys
# import os
import numpy as np
# import h5py
# import multiprocessing as mp
# import logging
import pdb

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3, venn2
import vcf
import re

def main():
    # read in the vcf files
    filename = 'T0diploid_S2_0_01.vcf'
    t = pd.read_table(filename, header = 1, skiprows = 9)
    T0POS = t.POS
    T0S1AF = map(AFfun, t.INFO)
    # S1AF = [float(s) for s in infonum]

    ## T1S3 diploid case; T0S2 diploid control
    filename = 'T1diploid_S3.vcf'
    T1T0POS, AF = table_import(filename)
    T0S2AF = AF[:,0]
    T1S3AF = AF[:,1]

    filename = 'T2diploid_S4.vcf'
    T2T0POS, AF = table_import(filename)
    T0S2AF_2 = AF[:,0]
    T2S4AF = AF[:,1]


    filename = 'T1diploid_S3_T2diploid_S4.vcf'
    T2T1POS, AF = table_import(filename)
    S3AF = AF[:,0]
    S4AF = AF[:,1]

    filename = 'T0haploid_S1.vcf'
    S2S1POS, AF = table_import(filename)
    S2AF = AF[:,0]
    S1AF = AF[:,1]


    ## Plot the positions as Venn diagram
    plt.figure(1)
    plt.subplot(212)
    # venn3([set(T1T0POS), set(T2T0POS),set(T2T1POS)], ('caseT1-controlT0_S2','caseT2-controlT0_S2','caseT2-controlT1'))
    venn3([set(T1T0POS), set(T2T0POS),set(T2T1POS)], ('T1-T0_S2','T2-T0_S2','T2-T1'))
    plt.subplot(211)
    # venn3([set(T0POS), set(T1T0POS), set(T2T0POS)], ('T0_S2','caseT1-controlT0_S2','T2-controlT0_S2'))
    venn3([set(T0POS), set(T1T0POS), set(T2T0POS)], ('T0_S2','T1-T0_S2','T2-T0_S2'))
    # venn3([set(T0POS), set(T1T0POS), set(T2T0POS),set(T2T1POS)], ('T0','T1-T0','T2-T0','T2-T1'))    
    plt.savefig('vennfig.png')
    # plt.figure(2)
    # y = [1,2,3]
    # for i in xrange(L):
    #     x = [T0AF[i]+i, T1AF[i]+i, T2AF[i]+i]
    #     plt.plot(x,y)
    # for i in xrange(L):
    #     plt.plot([i,i,i],[1,2,3],'--b')

    ## Plot the common positions for (T0 diploid control, T1 and T2 case) as bar plot
    T0T1T3POS = list(set(T1T0POS) & set(T2T0POS))
    T0T1T3POS = sorted(T0T1T3POS)
    L = len(T0T1T3POS)

    T0AF = [T0S2AF[list(T1T0POS).index(T0T1T3POS[i])] for i in xrange(L)] 
    T1AF = [T1S3AF[list(T1T0POS).index(T0T1T3POS[i])] for i in xrange(L)]
    T2AF = [T2S4AF[list(T2T0POS).index(T0T1T3POS[i])] for i in xrange(L)]

    ## barplot
    fig = plt.figure(2,figsize=(6,4))
    ax = fig.add_subplot(111)

    ## necessary variables
    ind = np.arange(L)                # the x locations for the groups
    width = 0.5                      # the width of the bars

    ## the bars
    rects1 = ax.bar(ind, T0AF, width/2.0,
                    color='black')
    rects2 = ax.bar(ind+width/2.0, T1AF, width/2.0,
                        color='red')
    rects3 = ax.bar(ind+width, T2AF, width/2.0,
                        color='blue')
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylabel('MAF (%)')
    ax.set_xlabel('Positions ChrXV780000+X')
    xTickMarks = [str(x-780000) for x in T0T1T3POS]
    ax.set_xticks(ind+width)
    xtickNames = ax.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=25, fontsize=10)
    ax.legend( (rects1, rects2, rects3), ('T0S2', 'T1','T2') ,prop={'size':10})
    plt.tight_layout()
    plt.savefig('T0S2T1T3POSbar.png')



    fig = plt.figure(3,figsize=(6,4))
    ax1 = fig.add_subplot(111)
    L = len(T1T0POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax1.bar(ind, T0S2AF, width,
                    color='black')
    rects2 = ax1.bar(ind+width, T1S3AF, width,
                        color='red')
    ax1.set_xlim(-width,len(ind)+width)
    ax1.set_ylabel('MAF (%)')
    xTickMarks = [str(x-780000) for x in T1T0POS]
    ax1.set_xticks(ind+width)
    xtickNames = ax1.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax1.legend( (rects1, rects2), ('control T0S2', 'case T1') ,loc=2,prop={'size':10})
    ax1.set_xlabel('Positions ChrXV:780000+X')
    plt.tight_layout()
    plt.savefig('T1_T0S2.png')


    fig = plt.figure(4,figsize=(6,4))
    ax2 = fig.add_subplot(111)
    L = len(T2T0POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax2.bar(ind, T0S2AF_2, width,
                    color='black')
    rects2 = ax2.bar(ind+width, T2S4AF, width,
                        color='red')
    ax2.set_xlim(-width,len(ind)+width)
    ax2.set_ylabel('MAF (%)')
    xTickMarks = [str(x-780000) for x in T2T0POS]
    ax2.set_xticks(ind+width)
    xtickNames = ax2.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax2.legend( (rects1, rects2), ('control T0S2', 'case T2'),loc=2,prop={'size':10} )
    ax2.set_xlabel('Positions ChrXV:780000+X')
    plt.tight_layout()
    plt.savefig('T2_T0S2.png')


    fig = plt.figure(5,figsize=(6,4))
    ax3 = fig.add_subplot(111)  
    L = len(T2T1POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax3.bar(ind, S3AF, width,
                    color='black')
    rects2 = ax3.bar(ind+width, S4AF, width,
                        color='red')
    ax3.set_xlim(-width,len(ind)+width)
    ax3.set_ylabel('MAF (%)')
    ax3.set_xlabel('Positions ChrXV:780000+X')
    xTickMarks = [str(x-780000) for x in T2T1POS]
    ax3.set_xticks(ind+width)
    xtickNames = ax3.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=25, fontsize=10)
    ax3.legend( (rects1, rects2), ('control T1', 'case T2'),loc=2 ,prop={'size':10})
    plt.tight_layout()    
    plt.savefig('T2_T1.png')


    fig = plt.figure(6,figsize=(6,4))
    ax3 = fig.add_subplot(111)  
    L = len(S2S1POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax3.bar(ind, S2AF, width,
                    color='black')
    rects2 = ax3.bar(ind+width, S1AF, width,
                        color='red')
    ax3.set_xlim(-width,len(ind)+width)
    ax3.set_ylabel('MAF (%)')
    ax3.set_xlabel('Positions ChrXV:780000+X')
    xTickMarks = [str(x-780000) for x in S2S1POS]
    ax3.set_xticks(ind+width)
    xtickNames = ax3.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=25, fontsize=10)
    ax3.legend( (rects1, rects2), ('control T0S2', 'case T0S1'),loc=2 ,prop={'size':10})
    plt.tight_layout()    
    plt.savefig('T0S2_T0S1.png')

    # pdb.set_trace()
    ## plot the three sets of experiments: 
    # T0 diploid control, T1 case; T0 diploid control, T2 case; T1 control, T2 case
    f, (ax1, ax2, ax3) = plt.subplots(3,figsize=(7,12))

    L = len(T1T0POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax1.bar(ind, T0S2AF, width,
                    color='black')
    rects2 = ax1.bar(ind+width, T1S3AF, width,
                        color='red')
    ax1.set_xlim(-width,len(ind)+width)
    ax1.set_ylabel('MAF (%)')
    xTickMarks = [str(x-780000) for x in T1T0POS]
    ax1.set_xticks(ind+width)
    xtickNames = ax1.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=25, fontsize=10)
    ax1.legend( (rects1, rects2), ('T0S2', 'T1') ,loc=2,prop={'size':10})


    L = len(T2T0POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax2.bar(ind, T0S2AF_2, width,
                    color='black')
    rects2 = ax2.bar(ind+width, T2S4AF, width,
                        color='red')
    ax2.set_xlim(-width,len(ind)+width)
    ax2.set_ylabel('MAF (%)')
    xTickMarks = [str(x-780000) for x in T2T0POS]
    ax2.set_xticks(ind+width)
    xtickNames = ax2.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=45, fontsize=10)
    ax2.legend( (rects1, rects2), ('T0S2', 'T2'),loc=2,prop={'size':10} )

    L = len(T2T1POS)
    ind = np.arange(L)
    width = 0.3
    rects1 = ax3.bar(ind, S3AF, width,
                    color='black')
    rects2 = ax3.bar(ind+width, S4AF, width,
                        color='red')
    ax3.set_xlim(-width,len(ind)+width)
    ax3.set_ylabel('MAF (%)')
    ax3.set_xlabel('Positions ChrXV780000+X')
    xTickMarks = [str(x-780000) for x in T2T1POS]
    ax3.set_xticks(ind+width)
    xtickNames = ax3.set_xticklabels(xTickMarks)
    plt.setp(xtickNames, rotation=25, fontsize=10)
    ax.legend( (rects1, rects2), ('T1', 'T2') )
    ax1.set_title('T0 diploid control, T1 case')
    ax2.set_title('T0 diploid control, T2 case')
    ax3.set_title('T1 control, T2 case')
    ax3.legend( (rects1, rects2), ('T1', 'T2'),loc=2 ,prop={'size':10})
    plt.tight_layout()

    plt.savefig('T0S2T1_T0S2T3_T1T3POSbar.png')

    # plt.show()

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