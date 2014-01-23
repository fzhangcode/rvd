
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
import rvd30

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)

    folderList = ('10','100','1000','10000')

    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        for f in folderList:
            path='vcf/%s' % str(100*(folderList.index(f))+1)
            if not os.path.exists(path):
                os.makedirs(path)
            controlFile = "./%s/Control.hdf5" %f
            caseFile = "Case%s.hdf5" % str(d).replace(".","_")
            caseFile = "./%(folder)s/%(file)s" %{'folder':f,'file':caseFile}
            print caseFile
            outputFile='%(path)s/vcf%(dilution)s' %{'path':path,'dilution':str(d).replace('.','_')}
           
            rvd30.test(controlFile,caseFile,T=0,N=1000,outputFile=outputFile)
if __name__ == '__main__':
    main()


"""import sys
import os
import rvd30
#Insert the source directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)

ControlFile = './Control.hdf5'
CaseFile = './Case0_1.hdf5'
OutputFile = 'vcf0_1'

assert os.path.isfile(ControlFile), "ControlFile file does not exist: %s" % ControlFile
assert os.path.isfile(CaseFile), "CaseFile file does not exist: %s" % CaseFile

rvd30.test(ControlFile, CaseFile, T=0, N=1000, outputFile=OutputFile)
"""
