import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import itertools as it

rvddir = os.path.join('./')
sys.path.insert(0, rvddir)
import rvd29
import rvd30
import rvd2

# Sensitivity analysis: Compare the M value in the depth chart = 10, dilution = 100%

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    N=1000 # Z sampling size  
    (n_gibbs, nmh) = (4000, 50)

# Plot E[M|Data] for RVD3 with Jeffreys prior
    logging.debug("Processing M for Jeffreys prior")
    caseFile = "./Case100_Jeffreys.hdf5" 
    (casePhi, caseM, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd29.load_model(caseFile)
    #fig = plt.figure()
    fig=plt.figure(figsize=(20, 12))
    ax1 = plt.subplot(131)
    caseLoc = [int(x.split(':')[1]) for x in caseLoc]
    #ax1.plot(caseLoc,caseN.T)
    ax1.semilogy(caseLoc,np.mean(caseM,axis=1),color='b',label="M")
    ax1.set_title('Jeffrey\'s prior')
    ax1.semilogy(caseLoc, list(it.repeat(casePhi['M0'],400)),linewidth=2,color='r',ls='--',label="M0")
    ax1.legend(loc='upper left')
    ax1.set_xlabel('Position',fontsize=15)
    ax1.set_ylabel('E(M|x)',fontsize=15)
    ax1.set_ylim(0.1, 100000000)
    ax1.plot([85,105,125,145,165,185,205,225,245,265,285,305,325,345],[1,1,1,1,1,1,1,1,1,1,1,1,1,1],'kx')


# Plot E[M|Data] for RVD3 with log-normal prior
    logging.debug("Processing M for log-normal prior")
    caseFile = "./Case100_lognormal.hdf5" 
    (casePhi, caseM, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd30.load_model(caseFile)
    ax2 = plt.subplot(132)
    caseLoc = [int(x.split(':')[1]) for x in caseLoc]
    #ax2.plot(caseLoc,caseN.T)
    ax2.semilogy(caseLoc,np.mean(caseM,axis=1),color='b')
    ax2.set_title('Log-normal prior',fontsize=15)
    ax2.semilogy([caseLoc[0],caseLoc[-1]],[casePhi['M0'],casePhi['M0']],linewidth=2,color='r',ls='--')
    ax2.set_xlabel('Position',fontsize=15)
    ax2.set_ylabel('E(M|x)',fontsize=15)
    ax2.set_ylim(0.1, 100000000)
    #ax2.set_yticks([1])
    #ax2.set_yticklabels(['M0'])
	
# Plot E[M|Data] for RVD3 with improper prior
    logging.debug("Processing M for improper prior")    
    caseFile = "./Case100_improper.hdf5" 
    (casePhi, caseTheta, caseMu, caseLoc, caseR, caseN) = rvd2.load_model(caseFile)
    ax3 = plt.subplot(133)
    caseLoc = [int(x.split(':')[1]) for x in caseLoc]
    #ax3.plot(caseLoc,caseN.T)
    ax3.semilogy(caseLoc,casePhi['M'],color='b')
    ax3.set_title('Improper prior')
    ax3.semilogy([caseLoc[0],caseLoc[-1]],[casePhi['M0'],casePhi['M0']],linewidth=2,color='r',ls='--')
    ax3.set_xlabel('Position',fontsize=15)
    ax3.set_ylabel('E(M|x)',fontsize=15)
    ax3.set_ylim(0.1, 100000000)
    #ax3.set_yticks([1])
    #ax3.set_yticklabels(['M0'])

    #plt.show()
    plt.savefig('priors_compare.pdf')

if __name__ == '__main__':
    main()

