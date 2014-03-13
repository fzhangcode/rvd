import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss

# Insert the src/python/rvd27 directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
import rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():

    controlFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Control.hdf5'
    caseFile = '../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome/Case.hdf5'
    
    hdf5name = 'HCC1187_call_unmasked.hdf5'

    try:
        with h5py.File(hdf5name, 'r') as f:
                somatic_call = f['somatic_call'][...]
                control_case_diff_call = f['control_case_diff_call'][...]
                case_control_diff_call =f['case_control_diff_call'][...]
    
                controlMu = f['controlMu'][...]
                caseMu = f['caseMu'][...]
                caseN_median = f['caseN_median'][...]
                controlN_median = f['controlN_median'][...]
                loc = f['loc'][...]
                f.close()

    except IOError as e:       
        diffalpha = 0.95
        diffroi = (0.0,np.inf)
        N = 1000
        outputFile = 'control_case_diff.vcf'
        [loc, control_case_diff_call, controlMu, caseMu, controlN, caseN, postP, chi2P ]= \
        rvd27.diff_test(diffalpha, controlFile, caseFile, diffroi, N, outputFile)
        
        ## median read depth
        caseN_median = np.median(caseN)
        controlN_median = np.median(controlN)
        

        outputFile = 'case_control_diff.vcf'
        [_,case_control_diff_call,_,_,_,_,_,_] = rvd27.diff_test(diffalpha, caseFile, controlFile, diffroi, N, outputFile)
        

        somatic_tau = (0.05, 0.75)
        [_,somatic_call,_,_,_,_] = rvd27.somatic_test(controlFile, caseFile, somatic_tau, 'somaticcall_all.vcf')
      
        with h5py.File(hdf5name, 'w') as f:
                f.create_dataset('somatic_call', data=somatic_call)
                f.create_dataset('control_case_diff_call', data=control_case_diff_call)
                f.create_dataset('case_control_diff_call', data=case_control_diff_call)
                f.create_dataset('caseN_median',data=caseN_median)
                f.create_dataset('controlN_median',data=controlN_median)
                f.create_dataset('caseMu',data=caseMu)
                f.create_dataset('controlMu',data=controlMu)
                f.create_dataset('loc',data=loc)
                f.close()
   
    call = [somatic_call,control_case_diff_call, case_control_diff_call]
    call = np.sum(call,0) > 0
    MuBarPlot(controlMu[call], caseMu[call], loc[call],'HCC1187_mu.pdf')

    print('Median read depth for control sample HCC1187 is %d' % controlN_median)
    print('Median read depth for case sample HCC1187 is %d' % caseN_median)
    

def MuBarPlot(CallControlMu,CallCaseMu,Loc,title):
    # 95% credible interval
    alpha = 0.05
    pos = int(CallControlMu.shape[1]*alpha/2)
    # pdb.set_trace()
    # sort along Gibbs samples
    ControlMu = np.sort(CallControlMu,axis=1)
    Control_yerr = np.array([np.mean(CallControlMu,1)-ControlMu[:,pos], ControlMu[:,CallControlMu.shape[1]-pos]-np.mean(CallControlMu,1)])

    CaseMu = np.sort(CallCaseMu,axis=1)
    Case_yerr = np.array([np.mean(CallCaseMu,1)-CaseMu[:,pos], CaseMu[:,CallCaseMu.shape[1]-pos]-np.mean(CallCaseMu,1)])
    
    if len(Loc) is not 0:
        ind = np.arange(CallCaseMu.shape[0])
        width=0.35
        fig, ax = plt.subplots(figsize=(9,6))

        rects1 = ax.bar(ind, np.mean(CallControlMu,1), width, color='0.8', yerr=Control_yerr,ecolor='k')
        rects2 = ax.bar(ind+width, np.mean(CallCaseMu,1), width, color='0.5', yerr=Case_yerr,ecolor='k')
        ax.set_xlim([0, 13])
    
        ax.set_ylabel('Minor Allele Frequency')
        ax.set_xlabel("HG19 Genomic Location [chr7:154,700,000+X]")
        ax.set_xticks(ind+width)
        
        diffcall = [1,2,3,5,6,8,9,10]
        for i in diffcall:
            ax.text(i+width-width/9.0, 1.0021,'*',color='r',fontsize=15)

        Loclabel = [x.split(':')[1][4:] for x in Loc]
        rslabel = ['rs1239326', ' ' ,'rs1239324', '       rs71534174', ' ' ,'rs35505514', ' ' , ' ' , ' ' , ' ' , ' ' ,'rs4398858', ' ']
 

        label = ['%s\n%s' %(Loclabel[i], rslabel[i]) for i in xrange(len(Loclabel))]

        ax.set_xticklabels(label)
        

        lgd = ax.legend( (rects1[0], rects2[0]), ('Control', 'Case'), bbox_to_anchor=(0.85, 1.02, 0.15, .102),
                   loc=4,ncol=1, mode="expand", borderaxespad=0., prop={'size':9})

        plt.rcParams.update({'font.size': 9, 'font.family': 'serif'})
        fig.tight_layout(rect=[0, 0, 0.95, 0.9])

        plt.savefig(title, bbox = 'tight', bbox_extra_artists=(lgd,))
        
    else:
        print 'No variant is called'

if __name__ == '__main__':
    main()
