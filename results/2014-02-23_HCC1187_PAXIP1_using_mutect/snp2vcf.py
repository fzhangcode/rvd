import sys
import os
import numpy as np
import logging
import pdb

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

        controlFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Control.hdf5'
        caseFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Case.hdf5'

        loc, refb, altb, controlMu, controlMu0, controlR, controlN, \
        caseMu, caseMu0, caseR, caseN = rvd27.load_dualmodel(controlFile, caseFile)

        position = np.array(['chr7:154743899', 'chr7:154749704', 
        'chr7:154753635', 'chr7:154754371', 'chr7:154758813',
         'chr7:154760439', 'chr7:154766700', 'chr7:154766732', 
         'chr7:154777014', 'chr7:154777118', 'chr7:154780960'])
        
        J = len(loc)
        call = []
        altb = []
        somatictype = []

        i = 0
        for j in xrange(J):
                if len(np.nonzero(loc[j] == position)[0]) ==0:
                        call.append(False)
                        altb.append('.')
                        somatictype.append('.')
                else:
                        call.append(True)
                        altb.append('.')
                        somatictype.append('.')
                        i+=1
        call = np.array(call)
        altb = np.array(altb)
        somatictype = np.array(somatictype)
        outputFile = 'HCC1187.vcf'
        somatic_tau = (0.05, 0.75)
        rvd27.write_dualvcf(outputFile,loc, call, refb, np.mean(controlMu, axis=1), \
              np.median(controlR,0), controlN, \
              np.mean(caseMu, axis=1), np.median(caseR,0), caseN, [somatic_tau], altb = altb, tag = 'somatic', somatictype = somatictype)
if __name__ == '__main__':
    main()
