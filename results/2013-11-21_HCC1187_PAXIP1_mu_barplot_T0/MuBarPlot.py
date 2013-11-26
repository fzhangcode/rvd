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
import rvd271 as rvd27

logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')

def main():
    #model samples
    
    controlFile = "../2013-11-07_HCC1187_PAXIP1_genome_fa/Control.hdf5"
    caseFile = "../2013-11-07_HCC1187_PAXIP1_genome_fa/Case.hdf5"
    
##    controlFile = "../2013-11-07_HCC1187_PAXIP1_hg19masked/Control.hdf5"
##    caseFile = "../2013-11-07_HCC1187_PAXIP1_hg19masked/Case.hdf5"
    
    hdf5name = 'HCC1187_call_unmasked.hdf5'
    
    try:
        with h5py.File(hdf5name, 'r') as f:
                call01 = f['call01'][...]
                call00 = f['call00'][...]
                call1 = f['call1'][...]
                call2 = f['call2'][...]
                call3 = f['call3'][...]
                call4 = f['call4'][...]
                call5 = f['call5'][...]
                call6 = f['call6'][...]
                call7 = f['call7'][...]
                call8 = f['call8'][...]
                controlMu = f['controlMu'][...]
                caseMu = f['caseMu'][...]
                caseN_median = f['caseN_median'][...]
                controlN_median = f['controlN_median'][...]
                loc = f['loc'][...]
                f.close()
    except IOError as e:
        
        diffalpha = 0.95
        diffroi = (0,np.inf)
        N = 1000
        outputFile = 'control_case_diff_plus.vcf'
        [loc, call00, controlMu, caseMu, controlN, caseN, postP, chi2P ]= rvd27.diff_test(0.95, controlFile, caseFile, \
                                                                     diffroi, N, outputFile = 'control_case_diff_plus.vcf')
        ## median read depth
        caseN_median = np.median(caseN)
        controlN_median = np.median(controlN)
        
        diffroi = (-np.inf, 0)
        [_,call01,_,_,_,_,_,_] = rvd27.diff_test(0.95, controlFile, caseFile, \
                                                diffroi, N, outputFile = 'control_case_diff_minus.vcf')
        Interval=[[0.0,0.2],[0.2,0.8],[0.8,1.0]]
        alpha = 0.8

        [_,call1,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[0],caseroi=Interval[1],outputFile='normal_case_ref_heter.vcf')
    ##    pdb.set_trace()
        [_,call2,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[0],caseroi=Interval[2],outputFile='normal_case_ref_homo.vcf')
        [_,call3,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[1],caseroi=Interval[0],outputFile='normal_case_heter_LOH1.vcf')
        [_,call4,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[1],caseroi=Interval[1],outputFile='normal_case_heter_unknown.vcf')
        [_,call5,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[1],caseroi=Interval[2],outputFile='normal_case_heter_LOH2.vcf')
        [_,call6,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[2],caseroi=Interval[0],outputFile='normal_case_homo_homo.vcf')
        [_,call7,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[2],caseroi=Interval[1],outputFile='normal_case_homo_heter.vcf')
        [_,call8,_,_,_,_,_,_] = rvd27.somatic_test(alpha, controlFile, caseFile, controlroi=Interval[2],caseroi=Interval[2],outputFile='normal_case_homo_unknown.vcf')    
        
        with h5py.File(hdf5name, 'w') as f:
                f.create_dataset('call01', data=call01)
                f.create_dataset('call00', data=call00)
                f.create_dataset('call1', data=call1)
                f.create_dataset('call2', data=call2)
                f.create_dataset('call3', data=call3)
                f.create_dataset('call4', data=call4)
                f.create_dataset('call5', data=call5)
                f.create_dataset('call6', data=call6)
                f.create_dataset('call7', data=call7)
                f.create_dataset('call8', data=call8)
                f.create_dataset('caseN_median',data=caseN_median)
                f.create_dataset('controlN_median',data=controlN_median)
                f.create_dataset('caseMu',data=caseMu)
                f.create_dataset('controlMu',data=controlMu)
                f.create_dataset('loc',data=loc)
                f.close()
       
    call = [call00,call01, call1,call2,call3,call4,call5,call6,call7,call8]
    print('Median read depth for control sample HCC1187 is %d' % controlN_median)
    print('Median read depth for case sample HCC1187 is %d' % caseN_median)
     
##    pdb.set_trace()
    call = np.sum(call,0) > 0
##    call = call00
    MuBarPlot(controlMu[call], caseMu[call], loc[call],'HCC1187_mu.pdf')

def MuBarPlot(CallControlMu,CallCaseMu,Loc,title):
    # 95% credible interval
    alpha = 0.05
    pos = int(CallControlMu.shape[1]*alpha/2)
    #pdb.set_trace()
    # sort along Gibbs samples
    ControlMu = np.sort(CallControlMu,axis=1)
    Control_yerr = np.array([np.mean(CallControlMu,1)-ControlMu[:,pos], ControlMu[:,CallControlMu.shape[1]-pos]-np.mean(CallControlMu,1)])

    CaseMu = np.sort(CallCaseMu,axis=1)
    Case_yerr = np.array([np.mean(CallCaseMu,1)-CaseMu[:,pos], CaseMu[:,CallCaseMu.shape[1]-pos]-np.mean(CallCaseMu,1)])
    
    if len(Loc) is not 0:
        ind = np.arange(CallCaseMu.shape[0])
        width=0.35
        fig, ax = plt.subplots()
##        rects1 = ax.bar(ind, np.mean(CallControlMu,1), width, color='r', yerr=Control_yerr,ecolor='k')
##        rects2 = ax.bar(ind+width, np.mean(CallCaseMu,1), width, color='g', yerr=Case_yerr,ecolor='k')
        rects1 = ax.bar(ind, np.mean(CallControlMu,1), width, color='0.8', yerr=Control_yerr,ecolor='k')
        rects2 = ax.bar(ind+width, np.mean(CallCaseMu,1), width, color='0.5', yerr=Case_yerr,ecolor='k')
        
        ax.set_ylabel('Minor Allele Frequency')
        ax.set_xlabel('Called Locations')
        ax.set_xticks(ind+width)

        ax.set_xticklabels( [x.split(':')[1] for x in Loc], rotation = 15)
##        ax.set_xticklabels( Loc, rotation = 60 )
        ax.set_ylim
        lgd = ax.legend( (rects1[0], rects2[0]), ('Control', 'Case'), bbox_to_anchor=(0., 1.02, 1., .102),
                   loc=3,ncol=2, mode="expand", borderaxespad=0. )
        plt.rcParams.update({'font.size': 15, 'font.family': 'serif'})
        plt.savefig(title, bbox = 'tight', bbox_extra_artists=(lgd,))
    else:
        print 'No variant is called'

if __name__ == '__main__':
    main()
