from __future__ import print_function
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp
import glob

import logging
import pdb
from datetime import date

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

def main():
    filelist = glob.glob('./output/10/calltbl*sample*.txt')

    foutput = 'stat.txt'
    statf = open(foutput,'w')
    print("filename\tTPR\tTNR\tFPR", file=statf)

    for f in filelist:
        [TPR, TNR, FDR ] = txt2vcf_single(f)
        print("%s\t%0.3f\t%0.3f\t%0.3f" %(f.split('/')[-1],TPR, TNR, FDR), file=statf)
    statf.close()

def txt2vcf_single(filename):
    callstatus= pd.read_table(filename, header=0)

    loc = np.array(callstatus.AlignReferencePosition)
    call = np.array(callstatus.Call)
    ref = np.arange(85,346,20)

    TP = 0
    TN = 0
    FP = 0
    FN = 0

    for i in xrange(len(loc)):
        t1 = call[i]
        t2 = np.nonzero(ref==loc[i])[0].size

        if t1 and t2:
            TP += 1
        elif not(t1) and not(t2): 
            TN += 1
        elif not(t1) and t2:
            FP += 1
        else:
            FP += 1

    if TP+FN != 0:
        TPR = float(TP)/(TP+FN)
    else:
        TPR=np.nan 
    #Specificity (TNR, true negative rate)
    if FP+TN != 0:
        TNR = float(TN)/(FP+TN)
    else:
        TNR = np.nan

    if FP+TP!=0:
        FDR = float(FP)/(FP+TP)
    else:
        FDR = np.nan



    return TPR, TNR, FDR 

if __name__ == '__main__':
    main()


