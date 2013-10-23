from __future__ import print_function
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import multiprocessing as mp

import logging
import pdb
from datetime import date

import vcf
# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    fracList = (10, 100, 1000, 10000)

    for f in fracList:
##        logging.debug("Processing dilution: %s" % str(d))
    	for d in dilutionList:
            logging.debug("Processing fraction: %s" % str(f))
            txtfile='work/%(frac)s/vcf%(dilution)s.txt' \
            %{'frac':str(f),'dilution':str(d).replace('.','_')}
##            pdb.set_trace()
            # import the call_stats file
##            callstatus= pd.read_table(txtfile, header=1, error_bad_lines=False)
            callstatus= pd.read_table(txtfile, header=31)
            chrom='1_1_400'
            pos=callstatus.POS
            ref=callstatus.REF
            alt=callstatus.ALT
            
            # write vcf file
            vcfdir='vcf/%s' % str(f)
            if not os.path.exists(vcfdir):
                os.makedirs(vcfdir)
            outputFile='%(vcfdir)s/vcf%(dilution)s.vcf' %{'vcfdir':vcfdir,'dilution':str(d).replace('.','_')}
            
            today=date.today()
            vcfF = open(outputFile,'w')
            
            print("##fileformat=VCFv4.1", file=vcfF)            
            print("##fileDate=%0.4d%0.2d%0.2d" % (today.year, today.month, today.day), file=vcfF)            
            print("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">", file=vcfF)
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAF", file=vcfF)
            for i in xrange(len(pos)):
                    print("%s\t%d\t.\t%c\t%s\t.\tPASS\t%0.1f" % (chrom, pos[i], ref[i], alt[i], 0), file=vcfF)            
            vcfF.close()

if __name__ == '__main__':
    main()
