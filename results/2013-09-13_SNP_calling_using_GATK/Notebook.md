2013-09-10 SNP calling using GATK
==============================

Purpose
------------
Compare GATK rvd2 on variants calling.

Conclusions
-----------------
GATK is able to detect variants at mutation rate 100% accurately with sensitivity/specificity = 1.00/1.00. GATK is also able to detect 10% snps with sensitivity at around 0.50 and specificity at 1.00. GATK is not able to detect variants if the mutation rate is 1% or lower.

Background
----------------
Please refer to the following webpage for the instructions.

[Calling SNPs with GATK’s Unified Genotyper](http://ged.msu.edu/angus/tutorials-2013/snp_tutorial.html#calling-snps-with-gatk-s-unified-genotyper)

[UnifiedGenotyper documentation](http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html)

[Introduction to the processing of short read next generation sequencing data](http://darrenjw.wordpress.com/2010/11/28/introduction-to-the-processing-of-short-read-next-generation-sequencing-data/)

Materials and Equipment
------------------------------
SNP_GATK.sh

Experimental Protocol
---------------------------
Run the Makefile in `../../data/2013-08-06_Downsample_Read_Depth` with command 'make -j 10' first to generate the downsampled sorted bam files. Bam files are downsampled at rate 0.1, 0.01, 0.001, 0.0001.

Run the SNP_GATK.sh to do SNP calling using GATK.


Results
-----------
The program generates vcf files which saves variant calling results. They are available in the same directory.

Using the characer.py in `../2013-09-19_operating_characteristics` it is possible to generate a summarizing table showing the sensitivity/specificity for variant calling using GATK. A screenshot is shown as follows:

![](http://i.imgur.com/dhac6U2.png)

Figure 1. Sensitivity/specificity achieved using GATK for variant calling.

It can be seen that GATK is able to detect variants at mutation rate 100% accurately with sensitivity/specificity = 1.00/1.00. GATK is also able to detect 10% snps with sensitivity at around 0.50 and specificity at 1.00. GATK is not able to detect variants if the mutation rate is 1% or lower.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 09/30/13


Witnessed by: ________________________
