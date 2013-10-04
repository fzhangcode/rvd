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
import logging
import pdb
import xlwt
# <codecell>

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(module)s:%(message)s')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

book=xlwt.Workbook(encoding="utf-8")
sheet1=book.add_sheet("Coverage")

dfracList=(10,100,1000,10000)
DilutionList = (0.1, 0.3, 1.0, 10.0,100.0)
for i in xrange(4):
    sheet1.write(0, i+2, "%s"%str(dfracList[i]))

for i in xrange(6):
    sheet1.write(i+1,0,"Control")
sheet1.write(7,0,"Control Median")
                 
for i in xrange(5):
    for j in xrange(6):
        sheet1.write(7*i+j+8, 0, "%0.1f%%" %DilutionList[i])
    sheet1.write(7*i+14, 0, "%0.1f%% Median" %DilutionList[i])

##pool =None
pool = mp.Pool(processes=48)
tocfilename = "../../data/synthetic_toc_p2.txt"
toc = pd.read_table(tocfilename)

controlFileList = ["../2013-08-06_Downsample_Read_Depth/depth_chart/10/%s" % filename for filename in toc.Filename[toc.isRef=='Y']]
for name in controlFileList:
    basename=os.path.basename(name)
    basename="%s" %basename
    sheet1.write(controlFileList.index(name)+1, 1,basename)

for dilution in DilutionList:
    caseFileList = ["../2013-08-06_Downsample_Read_Depth/depth_chart/10/%s" % filename for filename in toc.Filename[toc.Dilution==dilution]]
    for name in caseFileList:
        basename=os.path.basename(name)
        basename="%s" %basename
        sheet1.write(7*DilutionList.index(dilution)+caseFileList.index(name)+8, 1,basename)
    
for dfrac in dfracList:
    datapath="../2013-08-14_Compute_ROC_Synthetic_avg%d" %dfrac
    # import the coverage control
    logging.debug("Processing control data.")
    h5FileName = "%s/Control.hdf5" %datapath
##    pdb.set_trace()

    (_, _, _, _, _, controlN) = rvd27.load_model(h5FileName)
    Nmedian=np.median(controlN,1)
    for i in xrange(6):
        sheet1.write(i+1,dfracList.index(dfrac)+2,"%d" %Nmedian[i])
    Nmedian=np.median(Nmedian)
    sheet1.write(7,dfracList.index(dfrac)+2,"%d" %Nmedian)     

    # import the coverage for cases
    count=1
    for dilution in np.unique(toc[toc.isRef=='N'].Dilution):
        logging.debug("Processing dilution: %0.1f" % dilution)
        
        h5FileName = "Case%s.hdf5" % str(dilution).replace(".", "_", 1)
        h5FileName = "%(path)s/%(file)s" %{"path":datapath,"file":h5FileName}

        (_, _, _, _, _, caseN) = rvd27.load_model(h5FileName)
        Nmedian=np.median(caseN,1)
        for j in xrange(6):
            sheet1.write(7*count+j+1, dfracList.index(dfrac)+2, "%d" %Nmedian[j])
        Nmedian=np.median(Nmedian)
        sheet1.write(7*count+7, dfracList.index(dfrac)+2, "%d" %Nmedian)
        count=count+1
book.save('Converage.xls') 

