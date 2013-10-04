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
            workpath='work/%s' % str(f)
            txtfile='%(path)s/%(dilution)s.call_stats.txt' \
            %{'path':workpath,'dilution':str(d).replace('.','_')}

            # import the call_stats file
            callstatus= pd.read_table(txtfile, header=1)
            chrom=callstatus.contig
            pos=callstatus.position
            ref=callstatus.ref_allele
            var=callstatus.alt_allele
            
            # write vcf file
            outputFile='%(path)s/vcf%(dilution)s.vcf' %{'path':workpath,'dilution':str(d).replace('.','_')}
            
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
