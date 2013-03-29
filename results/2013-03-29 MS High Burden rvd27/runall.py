# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp

import subprocess

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd26 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27 as rvd

def main():
    bamFileNameList = [ r'../../data/ms/2012-09-12 Aligned Data/lane1_NoIndex_L001_R1.fastq.gz.prefix.bam',
                       r'../../data/ms/2012-09-12 Aligned Data/lane2_NoIndex_L002_R1.fastq.gz.prefix.bam',
                       r'../../data/ms/2012-09-12 Aligned Data/lane3_NoIndex_L003_R1.fastq.gz.prefix.bam']
    # bamFileNameList = [ r'../../data/ms/2012-09-12 Aligned Data/lane5_NoIndex_L005_R1.fastq.gz.prefix.bam',
    #                     r'../../data/ms/2012-09-12 Aligned Data/lane6_NoIndex_L006_R1.fastq.gz.prefix.bam',
    #                     r'../../data/ms/2012-09-12 Aligned Data/lane7_NoIndex_L007_R1.fastq.gz.prefix.bam']

    logging.info("Converting BAM files to PILEUP.")
    pileupFileNameList = [ make_pileup(bamFileName,
                                 r'../../data/hg19/hg19.fa.masked',
                                 r'chr18:67529192-67625232') 
                           for bamFileName in bamFileNameList]
    
    logging.info("Converting PILEUP files to DEPTH CHART.")
    dcFileNameList = [ make_depth(pileupFileName)
                       for pileupFileName in pileupFileNameList]
    
    logging.info("Loading depth charts.")
    (r, n, loc, refb) = load_depth(dcFileNameList)
    
    # Apply a minimum depth threshold
    L=10
    depthInd = np.sum(n, axis=0) > L
    logging.debug("Number of positions before depth filter at L=%d: %d. After: %d." % ( L, np.shape(n)[1], np.sum(depthInd) ) )
    r = r[:,depthInd]
    n = n[:,depthInd]

    # pool = mp.Pool(processes=62)
    pool = None    

    h5FileName = "control.hdf5"
    try:
        with h5py.File(h5FileName, 'r') as f:
            pass
    except IOError as e:
        phi, theta_s, mu_s = rvd.mh_sample(r, n, nsample=1000, burnin=0.2, pool=pool)
        logging.debug("Saving model in %s" % h5FileName)
        rvd.save_model(h5FileName, loc, refb, r, n, phi, theta_s, mu_s)
    
def load_depth(dcFileNameList):
    """ Return (r,n) for a list of depth charts.
    """
    r=[]; n=[]
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    loc = []
    refb = {}
    cd = []
    for dcFileName in dcFileNameList:
        with open(dcFileName, 'r') as dcFile:
            header = dcFile.readline().strip()
            dc = dcFile.readlines()
            dc = [x.strip().split("\t") for x in dc]
    
            loc1 = map(int, [x[2] for x in dc if x[4] in acgt.keys()])
            loc.append( loc1 )
            
            refb1 = dict(zip(loc1, [x[4] for x in dc if x[4] in acgt.keys()]))
            refb.update(refb1)
            cd.append( dict(zip(loc1, [map(int, x[5:9]) for x in dc if x[4] in acgt.keys()])) )
            
    loc = list(reduce(set.intersection, map(set, loc)))
    refb = [refb[k] for k in loc]
    
    J = len(loc)
    N = len(dcFileNameList)
    c = np.zeros( (J, 4, N) )
    for i in xrange(0, N):
            c = np.array( [cd[i][k] for k in loc] )
            n1 = np.sum(c, 1)
            r1 = np.zeros(J)
            for j in xrange(0,J):
                r1[j] = n1[j] - c[j, acgt[refb[j]]]
            r.append(r1)
            n.append(n1)
    r = np.array(r)
    n = np.array(n)

    return (r,n,loc, refb)



def make_depth(pileupFileName):
    if not os.path.isdir("depth_chart"):
        os.makedirs("depth_chart")
        
    callString = ["../../bin/pileup2dc", "%s" % pileupFileName]

    dcFileName = pileupFileName.split("/")[-1]
    dcFileName = os.path.join("depth_chart", 
                                  "%s.dc" % dcFileName.split(".", 1)[0])
    try:
        with open(dcFileName, 'r'):
		logging.debug("Depth chart file exists: %s" % dcFileName) 
    except IOError:
        logging.debug("Converting %s to depth chart." % pileupFileName)
        with open(dcFileName, 'w') as fout:
            subprocess.call(callString, stdout=fout)
    return dcFileName
    
def make_pileup(bamFileName, fastaFileName, region):
    """ Creates a pileup file using samtools mpileup in /pileup directory.
    """
    
    # Check that BAM file exists
    assert os.path.isfile(bamFileName), "BAM file does not exist: %s" % bamFileName
    
    # Check that FASTA reference file exists
    assert os.path.isfile(fastaFileName), "FASTA file does not exist: %s" % fastaFileName
    
    # Create pileup directory if it doesn't exist
    if not os.path.isdir("pileup"):
        os.makedirs("pileup")
    
    # Format the samtools call
    callString = ["samtools", "mpileup", "-d", "1000000", "-r", "%s" % region,
                  "-f", "%s" % fastaFileName, "%s" % bamFileName]
                  
    # Remove the extension from the bam filename and replace with .pileup
    pileupFileName = bamFileName.split("/")[-1]
    pileupFileName = os.path.join("pileup", 
                                  "%s.pileup" % pileupFileName.split(".", 1)[0])
    
    # Run samtools pileup only if the file doesn't already exist.
    try:
        with open(pileupFileName, 'r'):
            logging.debug("Pileup file exists: %s" % pileupFileName)
    except IOError:
        logging.debug("[call] %s", " ".join(callString))
        with open(pileupFileName, 'w') as fout:
            subprocess.call(callString, stdout=fout)
    return pileupFileName
if __name__ == '__main__':
    main()
