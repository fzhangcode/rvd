from __future__ import print_function
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp
# <codecell>
import logging
import pdb
from datetime import date

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    fracList = (10, 100, 1000, 10000)

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
    	for f in fracList:
            vcfpath='vcf/%s' % str(f)
            if not os.path.exists(vcfpath):
                os.makedirs(vcfpath)

            snppath='snp/%s' % str(f)
            snpfile='%(path)s/snp%(dilution)s.snp' %{'path':snppath,'dilution':str(d).replace('.','_')}

            # import snp file
            chrom=[]
            pos=[]
            ss=[]
            ref=[]
            var=[]
            with open(snpfile,'rb') as source:
                for line in source:
                    fields = line.split('\t')
                    chrom.append(fields[0])
                    pos.append(fields[1])        
                    ref.append(fields[2])
                    var.append(fields[3])
                    ss.append(fields[12])

            idx=[i for i, x in enumerate(ss) if x=='Somatic' or x=='LOH' or x=='IndelFilter' or x=='Unknown']
            idx=np.array(idx)

            chrom=np.array([chrom[i] for i in idx])
            pos=np.array([int(pos[i]) for i in idx])
            ref=np.array([ref[i] for i in idx])
            var=np.array([var[i] for i in idx])
            
            # write vcf file
            outputFile='%(vcfpath)s/vcf%(dilution)s.vcf' %{'vcfpath':vcfpath,'dilution':str(d).replace('.','_')}
            
            today=date.today()
            vcfF = open(outputFile,'w')
            
            print("##fileformat=VCFv4.1", file=vcfF)            
            print("##fileDate=%0.4d%0.2d%0.2d" % (today.year, today.month, today.day), file=vcfF)            
            print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">", file=vcfF)
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAF", file=vcfF)
            for i in xrange(len(pos)):
                    print("%s\t%d\t.\t%c\t%s\t.\tPASS\t%0.1f" % (chrom[i], pos[i], ref[i], var[i], 0), file=vcfF)            
            vcfF.close()

if __name__ == '__main__':
    main()
