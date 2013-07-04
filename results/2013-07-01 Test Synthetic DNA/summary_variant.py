#!/usr/bin/env python
""" 
"""

import sys
import os
import numpy as np

# Insert the src/python/rvd26 directory at the front of the path.
base_dir = os.path.dirname(__file__)
rvddir = os.path.join(base_dir, '../../src/python/rvd26')
sys.path.insert(0, rvddir)
import rvd26

def main():
    filename = "../../data/synthetic_toc.txt"
    dcFileList = process_toc(filename)
    
    r=[]; n=[]
    
    for dcFile in dcFileList:
        dcFile = "../../data/synthetic_dcs/%s" % dcFile
        (r1, n1) = load_depth(dcFile)
        r.append(r1)
        n.append(n1)
    r = np.array(r)
    n = np.array(n)
    
    phi, theta, mu, M = rvd26.gibbs_map(r, n)
    

def load_depth(filename):
    acgt = {'A':0, 'C':1, 'G':2, 'T':3}
    
    dcFile = open(filename, 'r')
    header = dcFile.readline().strip()
    dc = dcFile.readlines()
    dc = [x.strip().split("\t") for x in dc]
    
    loc = map(int, [x[2] for x in dc])
    refb = [x[4] for x in dc]
    c = [map(int, x[5:9]) for x in dc]
    c = np.array(c)
    
    (J, K) = np.shape(c)
    n = np.sum(c, 1)
    r = np.zeros(J)
    for j in xrange(0, J):
        r[j] = n[j] - c[j, acgt[refb[j]]]
    return (r, n)

def process_toc(filename):
    tocfile = open(filename,'r')
    header = tocfile.readline().strip("\t")
    
    dcFile = []
    for i in xrange(0,3):
        (lane, pair, tagindex, tag, dilution, isRef, file) \
            = tocfile.readline().strip().split("\t")
        dcFile.append(file)
    
    return dcFile

if __name__ == '__main__':
    main()
    
    