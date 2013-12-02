
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the source directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd29

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('./')

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        for f in folderList:
            #path='vcf/%s' % str(10**(folderList.index(f)+1))
            path='vcf/%s' % str(10*(folderList.index(f)))
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "./%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "./%(file)s" %{'folder':f,'file':caseFile}
            
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
           
            rvd29.test(controlFile,caseFile,T=0,N=1000,outputFile=outputFile)
if __name__ == '__main__':
    main()


"""import sys
import os
import rvd29
#Insert the source directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)

ControlFile = './Control.hdf5'
CaseFile = './Case0_1.hdf5'
OutputFile = 'vcf0_1'

assert os.path.isfile(ControlFile), "ControlFile file does not exist: %s" % ControlFile
assert os.path.isfile(CaseFile), "CaseFile file does not exist: %s" % CaseFile

rvd29.test(ControlFile, CaseFile, T=0, N=1000, outputFile=OutputFile)
"""
