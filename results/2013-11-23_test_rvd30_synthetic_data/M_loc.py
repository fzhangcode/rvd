import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb


# Insert the rvd30 directory at front of the path
rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd30

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    dilutionList = (0.1,0.3,1.0,10.0,100.0)
    
    folder = '2013-11-23_test_rvd30_synthetic_data'
    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)

    fig=plt.figure(figsize=(12,20))
    #plt.suptitle('Read depth/M across position')

    controlFile = "../%s/Control.hdf5" %folder
    (controlPhi, controlM, controlTheta, controlMu, controlLoc, controlR, controlN) = rvd30.load_model(controlFile)

##  pdb.set_trace()
    sub0=len(dilutionList)+1
    ax=fig.add_subplot(sub0,2,1)
    #TODO: use index of controlN rather than directly controlLOC
    controlLoc = [int(x.split(':')[1]) for x in controlLoc]
    ax.plot(controlLoc,controlN.T)
    ax.set_title('Control')
    ax.set_ylabel('Coverage')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    ax=fig.add_subplot(sub0,2,2)
    ax.semilogy(controlLoc,np.mean(controlM,axis=1))
    ax.set_title('Control')
    ax.set_ylabel('M')
    ax.semilogy([controlLoc[0],controlLoc[-1]],[controlPhi['M0'],controlPhi['M0']],color='r',ls='--')
    for d in dilutionList:
        logging.debug("Processing dilution: %0.1f%%" % d)
        caseFile = "Case%s.hdf5" % str(d).replace(".","_")
        caseFile = "../%(folder)s/%(file)s" %{'folder':folder,'file':caseFile}
        (casePhi, caseM, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd30.load_model(caseFile)
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+3)
        caseLoc = [int(x.split(':')[1]) for x in caseLoc]
        ax.plot(caseLoc,caseN.T)
        ax.set_title("Dilution %0.1f%%" % d)
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Position')
        ax.set_ylabel('Coverage')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax=fig.add_subplot(sub0,2,2*dilutionList.index(d)+4)
        ax.semilogy(caseLoc,np.mean(caseM,axis=1))
        ax.set_title("Dilution %0.1f%%" % d)
        if dilutionList.index(d)==len(dilutionList)-1:
            ax.set_xlabel('Position')
        ax.set_ylabel('M')
##        pdb.set_trace()   
        ax.semilogy([caseLoc[0],caseLoc[-1]],[casePhi['M0'],casePhi['M0']],color='r',ls='--')
    plt.savefig('M_loc.pdf')

if __name__ == '__main__':
    main()

