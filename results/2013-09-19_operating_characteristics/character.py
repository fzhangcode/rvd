
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import logging
import pdb
import scipy.stats as ss
import vcf
import xlwt


def main():
    book=xlwt.Workbook(encoding="utf-8")
    sheet1=book.add_sheet("Sheet1")
    sheet1.write(0, 0, "Dilution")
    sheet1.write(0, 1, "Depth")

    sheet2=book.add_sheet("Sheet2")
    sheet2.write(0, 0, "Dilution")
    sheet2.write(0, 1, "Depth")
    
    method = {'rvd2':'../2013-08-15_Compute_ROC_Synthetic_avg_all/vcf', 'samtools':'../2013-09-10_SNP_calling_using_samtools/vcf','GATK':'../2013-09-13_SNP_calling_using_GATK/vcf'}
    
    DilutionList = (0.1, 0.3, 1.0, 10.0)
    DepthList = (10, 100, 1000, 10000)
    i=0
    
    for k, v in method.iteritems():
        i=i+1
        sheet1.write(0, i+1, k)
        sheet2.write(0, 9*(i-1)+6, k)
        character=('Sensitiviy', 'Specificity', 'FPR', 'FNR', 'PPV', 'NPV', 'FDR', 'ACC', 'MCC')
        for j in xrange(9):
            sheet2.write(1,9*(i-1)+j+2,character[j])
        for d in DilutionList:
            if i==1:
                sheet1.write(DilutionList.index(d)*len(DepthList)+1,0,"%0.1f%%" %d)
                sheet2.write(DilutionList.index(d)*len(DepthList)+2,0,"%0.1f%%" %d)
            for r in DepthList:
                if i==1:
                    sheet1.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1, 1, "%s" % str(r))
                    sheet2.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+2, 1, "%s" % str(r))
                    print DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1
                vcfFile=os.path.join(v,"%s" %r,
                                     "vcf%s.vcf" %str(d).replace('.','_'))
                vcf_reader = vcf.Reader(open(vcfFile, 'r'))
                callpos=np.array([record.POS for record in vcf_reader])
                
                RefClass = np.zeros(400)
                pos = np.arange(85,346,20)
                RefClass[pos-1] = np.ones_like(pos)

                PredictClass = np.zeros(400)
               # pdb.set_trace()
                if len(callpos) != 0:
                    PredictClass[callpos-1] = np.ones_like(callpos)
                
                [TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC]=characteristics(RefClass, PredictClass)
                ncharacter=(TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC)
                sheet1.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+1,i+1,"%(TPR)0.2f/%(TNR)0.2f" %{'TPR':TPR,'TNR':TNR})
                for j in xrange(9):
                    if np.isnan(ncharacter[j]):
                        sheet2.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+2,9*(i-1)+j+2,'NaN')
                    else:
                        sheet2.write(DilutionList.index(d)*len(DepthList)+DepthList.index(r)+2,9*(i-1)+j+2,ncharacter[j])
    book.save('character.xls')  

def characteristics(RefClass = None, PredictClass = None):

    if RefClass is None:
        RefClass = np.zeros(400)
        pos = np.arange(85,346,20)
        RefClass[pos-1] = np.ones_like(pos)
    if PredictClass is None:
        PredictClass = np.copy(RefClass)
        pos = np.arange(85,246,20)
        PredictClass[pos-1] = np.zeros_like(pos)



    #True Positive
    TP = len([i for i in range(len(RefClass)) if RefClass[i]==1 and PredictClass[i]==1])
    #True Negative
    TN = len([i for i in range(len(RefClass)) if RefClass[i]==0 and PredictClass[i]==0])
    #False Positive
    FP = len([i for i in range(len(RefClass)) if RefClass[i]==0 and PredictClass[i]==1])
    #False Negative
    FN = len([i for i in range(len(RefClass)) if RefClass[i]==1 and PredictClass[i]==0])
  
    #RefClassence Postive
    P1 = sum(RefClass)
    #RefClassence Negative
    N1 = len(RefClass) - P1

    #PredictClassion Postive
    P2 = sum(PredictClass)
    #PredictClassion Negative
    N2 = len(PredictClass) - P2

    #Sensitivity (TPR,true positive rate)
    if TP+FN != 0:
        TPR = float(TP)/(TP+FN)
    else:
        TPR=np.nan 
    #Specificity (TNR, true negative rate)
    if FP+TN != 0:
        TNR = float(TN)/(FP+TN)
    else:
        TNR = np.nan

    #FPR (1-Specificity, false negative rate)
    if FP+TN != 0:
        FPR = float(FP)/(FP+TN)
    else:
        FPR = np.nan
    #FNR
    if TP+FN !=0:
        FNR = float(FN)/(TP+FN)
    else:
        FNR = np.nan

    #PPV (Positive PredictClassive value)
    #pdb.set_trace()
    if TP+FP !=0:
        PPV = float(TP)/(TP+FP)
    else:
        PPV = np.nan
    #NPV (Negative PredictClassive value)
    if TN+FN !=0:
        NPV = float(TN)/(TN+FN)
    else:
        NPV = np.nan

    #FDR (false discovery rate)
    if FP+TP!=0:
        FDR = float(FP)/(FP+TP)
    else:
        FDR = np.nan

    #ACC (Accuracy)
    ACC = (TP+TN)/(P1+N1)
    #MCC (Matthews correlation coefficient)
    if P1*N1*P2*N2 != 0:
        MCC = float(TP*TN-FP*FN)/np.sqrt(P1*N1*P2*N2)
    else:
        MCC= np.nan
    # AUC

    return  TPR, TNR, FPR, FNR, PPV, NPV, FDR, ACC, MCC

if __name__ == '__main__':
    main()
