2013-09-19 Statistical measures table for the performance of various mutation detection method. 
=================================

Purpose
------------
To generate the statistical measures table for the performance of various mutation detection method.

Conclusions
-----------------
Result table is in the excel file character.xls in the current folder.
Background
----------------


Materials and Equipment
------------------------------
vcf files in the folders provided in Experimental Protocal section.

character.py in the current folder to generate various statistical characteristics.
  
Experimental Protocol
---------------------------
Method evaluated and corresponding vcf files:

	'rvd2_optimalT':'../2013-08-15_Compute_ROC_Synthetic_avg_all/vcf',
    'rvd2_half_dilution':'../2013-09-21_SNP_calling_RVD2_half_dilution/vcf',
    'rvd2_zero':'../2013-09-27_SNP_calling_RVD2_zero/vcf',
    'VarScan2_somatic':'../2013-09-23_SNP_calling_using_varscan2_somatic/vcf',
    'samtools':'../2013-09-10_SNP_calling_using_samtools/vcf',
    'GATK':'../2013-09-13_SNP_calling_using_GATK/vcf',
    'VarScan2':'../2013-09-20_SNP_calling_using_varscan2/vcf'.

**rvd2_optimalT**: statistical measures obtained from rvd2 approach with threshold set as optimal threshold returned from roc curves.

**rvd2_half_dilution**: statistical measures obtained from rvd2 approach with threshold set at half dilution rate.

**rvd2_zero**: statistical measures obtained from rvd2 approach with threshold set at zero.

**VarScan2**: VarScan2 version 2.3.2 is used via a pipe from samtools mpileup. In the experiment, -C50 option is specified because of the using of BWA form alignment. -d was set at 100000 to ensure the high read depth in our dataset. The minimum variant frequency was set to 0.00001 according to review paper [Accurately identifying low-allelic fraction variants in single samples with next-generation sequencing: applications in tumor subclone resolution](http://www.ncbi.nlm.nih.gov/pubmed/23766071)

The manual Page for VarScan2 mutation detection is in [mpileup2snp](http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp)

Please refer to the paper [VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing](http://www.ncbi.nlm.nih.gov/pubmed/22300766) for more information about the algorithm. 

**VarScan2_somatic**: The somatic mode of VarScan2 calls variants and identifies their somatic status. Please refer to manual page at [Varscan somatic](http://varscan.sourceforge.net/using-varscan.html#v2.3_somatic). For the parameters setting, normal-purity is set to be equal to dilution rate. min-var-freq is set at 0.00001. Other parameters are set as default.

**Samtools**: Please refer to web page 
[Calling SNPs/INDELs with SAMtools/BCFtools](http://samtools.sourceforge.net/mpileup.shtml )
and
[Calling SNPs with Samtools](http://ged.msu.edu/angus/tutorials-2013/snp_tutorial.html ) for more information about variant detection using Samtools.

**GATK**:Please refer to the following webpage for the instructions.[Calling SNPs with GATK’s Unified Genotyper](http://ged.msu.edu/angus/tutorials-2013/snp_tutorial.html#calling-snps-with-gatk-s-unified-genotyper) and [UnifiedGenotyper documentation](http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html).

Various measures are computed to evaluate the performance of these different methods.

So for a experiment with P positive instances and N negative instances, based on a certain classification method the four outcomes can be formulated in a 2*2 contigency table as:

![](http://i.imgur.com/YeJBq4S.png)

Different statistics are defined as:

![](http://i.imgur.com/ED9zXUY.png)

Please refer to web page [Receiver Operating Characteristics](http://en.wikipedia.org/wiki/Receiver_operating_characteristic) for more information.
Results
-----------
Result table is in the excel file character.xls in the current folder. The excel book has two sheets. The data cell in sheet1 shows sensitivity and specificity with the format of sensitivity/specificity. Sheet2 contains a variety of measures defined previously. A screenshot of sheet1 is shown as following:
![](http://i.imgur.com/3tcK2U4.png)

From the table it can be seen that

1. All the methods except rvd2 with threhold at half dilution rate were able to detect variants at 100% mutation rate. rvd2_optimal, rvd2_zero, VarScan2 and GATK achieved 1.00 sensitivity and 1.00 specificity when detecting 100% mutated variants.
 
2. Samtools, VarScan2 mpileup2snp were not able to detect variants when then mutation rate is 10%, 1.0%, 0.3% and 0.1% (Sensitivity=0). GATK showed higher power the samtools and VarScan2 mpileup2snp and achieved sensitivity=0.43,0.57,0.57 and 0.79, specificity=1 when the mutation rate is at 10.0%. But GATK also have no detection power regarding to 1.0%, 0.3%, 0.1% variants (sensitivity=0).

3. VarScan2_somatic has problem with high read depth, which is against the general rule of statistics. One thing worth noting is that for VarScan2 in somatic mode, the algorithm usually performed better when the read depth is low than when the read depth is high. VarScan2_somatic achieved sensitivity=1 and specificity=0 for 100% dilution rate when coverage was at 30590, a bad performance comparing to any other methods. Also, when the dilution rate is 10% and coverage at 26959, the sensitivity/specificity is 0.50/0.49.
 
4. Among all the methods, rvd2 with optimal threshold outcompeted other methods basically across all dilution rate and coverage. rvd2_optimalT achieved sensitivity=1.00 and specificity=1.00 when dilution rate is 10% and 100% across almost all coverage. The only exception is when dilution rate is 10% and coverage average is 22, the sensitivity is 0.86 and specificity is 1.00, which means rvd2 fails to detect 14% True positive. For dilution rate 1.0% and 0.3%, rvd2_optimalT needs coverage at 10 to the power of 3 to achieve decent sensitivity and specificity (higher than 0.90). For dilution rate at 0.1%, rvd2_optimal achieve sensitivity/specificity=0.71/0.76 when coverage average is 4129 and sensitivity/specificity=0.93/0.88 when coverage average is 41449. 

6. Since the optimalT is obtained from the ROC curves, which is not available for the real data, rvd2_zero set threshold at zero for all dilution rate and all read depth. rvd2_zero achieved relatively worse outcome comparing to optimalT, but still outcompete other approaches including VarScan2_somatic. for low dilution rate like 0.1% and 0.3%, rvd2 needs coverage up to 10 to power 4 to achieve decent sensitivity and specificity(higher or around 0.90), but rvd2_zero achieved sensitivity/specificity=1.00/1.00 for all read depth when dilution rate at 100% and rvd2_zero achieved sensitity/specificity=1.00/1.00 for dilution rate at 10% except when read depth is too low at 22. Coverage at 535 is enough for rvd2_zero to achieve sensitivity/specificity=0.93/0.95 when dilution rate is 1.0%.

7. rvd2 achieve sensitivity at zero when threshold set at half dilution. Therefore, half dilution is too high a threshold for all dilution rate. 

 
Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He_________________     Date: ________09/30/13_____________


Witnessed by: ________________________
