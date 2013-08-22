import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb

##os.chdir('S://yhe2//Research//rvd2//results//2013-08-21_depth_M_plot')

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    
    folder = '2013-08-19_Compute_ROC_Synthetic_avg10000'
    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)

    fig=plt.figure(figsize=(12,20))
    plt.suptitle('Read depth and M plot')

    controlFile = "../%s/Control.hdf5" %folder
    (controlPhi, controlTheta, controlMu, controlLoc, controlR, controlN) = rvd27.load_model(controlFile)
##    pdb.set_trace()
    
    sub0=len(dilutionList)+1
    ax=fig.add_subplot(sub0,2,1)
    ax.plot(controlLoc,controlN.T)
    ax.set_title('Control')
    ax.set_ylabel('Depth(N)')
##    leg=['replicate1','2','3','4','5','6']
##    ax.legend(leg)
    ax=fig.add_subplot(sub0,2,2)
    ax.plot(controlLoc,controlPhi['M'])
    ax.set_title('Control')
    ax.set_ylabel('M')
    
    
    
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
        caseFile = "../%(folder)s/%(file)s" %{'folder':folder,'file':caseFile}
        (casePhi, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd27.load_model(caseFile)
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+3)
        ax.plot(caseLoc,caseN.T)
        ax.set_title(str(d).replace(".","_"))
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Loc')
        ax.set_ylabel('Depth(N)')
        
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+4)
        ax.plot(caseLoc,casePhi['M'])
        ax.set_title(str(d).replace(".","_"))
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Loc')
        ax.set_ylabel('M')
    plt.savefig('Depth_M')
    plt.show()  

if __name__ == '__main__':
    main()

