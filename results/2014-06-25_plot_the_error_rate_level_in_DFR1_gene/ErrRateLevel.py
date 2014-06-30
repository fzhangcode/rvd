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

    ########################################### Plot the MAF level of all the positions#################################
    filename = './../2014-06-18_apply_rvd27_to_elisa_data/hdf5folder/T0diploid_S2.hdf5'
    MAFlevel_plot(filename, stitle = 'T0diploid_S2.png')
    
    filename = './../2014-06-18_apply_rvd27_to_elisa_data/hdf5folder/T0haploid_S1.hdf5'
    MAFlevel_plot(filename, stitle = 'T0haploid_S1.png')

    filename = './../2014-06-18_apply_rvd27_to_elisa_data/hdf5folder/T1diploid_S3.hdf5'
    MAFlevel_plot(filename, stitle = 'T1diploid_S3.png')

    filename = './../2014-06-18_apply_rvd27_to_elisa_data/hdf5folder/T2diploid_S4.hdf5'
    MAFlevel_plot(filename, stitle = 'T2diploid_S4.png')    

def MAFlevel_plot(filename, stitle = 'MAF.png'):
    with h5py.File(filename,'r') as f:
        mu = f['mu'][...]
        loc = f['loc'][...]
    loc = [int(x.split(':')[1]) for x in loc]

    mu = np.mean(mu, 1)
    L = len(loc)
    fig = plt.figure()
    ax = fig.add_subplot(111) 

    idx = np.arange(L) +1
    ax.plot(idx, mu*1000)

    ax.set_xlim(0, 650)
    ax.set_ylabel(r'$\mu [0/00]$')
    ax.set_xlabel('Position index')


    ## Label the highest ten positions
    sortidx = np.argsort(mu)

    for i in sortidx[-10:]:
        ax.text(idx[i], mu[i]*1000,str(loc[i]), rotation=90, va="bottom", ha="center",)
    plt.savefig(stitle)

if __name__ == "__main__":
    main()