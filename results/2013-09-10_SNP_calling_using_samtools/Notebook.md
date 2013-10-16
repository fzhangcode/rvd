2013-09-10 SNP calling using SAMtools/BCFtools
==============================

Purpose
------------
Using SAMtools mpileup and BCTtools to call variants.

Conclusions
-----------------
SAMtools is able to detect variants at mutation rate 100% and SAMtools is not able to detect variants if the mutation rate is 10% or lower. 

Background
----------------
Please refer to the following webpage for the instructions.

[Calling SNPs/INDELs with SAMtools/BCFtools](http://samtools.sourceforge.net/mpileup.shtml )
and
[Calling SNPs with Samtools](http://ged.msu.edu/angus/tutorials-2013/snp_tutorial.html )
Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
Run the Makefile in `../../data/2013-08-06_Downsample_Read_Depth` with command 'make -j 10' first to generate the downsampled sorted bam files. Bam files are downsampled at rate 0.1, 0.01, 0.001, 0.0001.

Run the SNP_samtools.sh to do SNP calling using samtools/bcftools.
samtools mpileup is used for variant calling. -C50 option is specified because of the using of BWA to form alignment. -d was set at 100000 to ensure the high read depth in our dataset. For each dilution rate and each downsampling rate, 6 bam files are fed into samtools mpile up as input. The reference fasta file is from `../../data/Synthetic_BAM_files/plasmid.fa`


Results
-----------
The program generates vcf files which saves variant calling results. They are available in the same directory.

Using the characer.py in `../2013-09-19_operating_characteristics` it is possible to generate a summarizing table showing the sensitivity/specificity for variant calling using SAMtools. A screenshot is shown as follows:

![](http://i.imgur.com/SbfEVFt.png)

Figure 1. Sensitivity/specificity achieved using samtools for variant calling.

It can be seen that SAMtools is able to detect variants at mutation rate 100% and SAMtools is not able to detect variants if the mutation rate is 10% or lower. Samtools can achieve sensitivity/specificity = 1.00/0.99 when the read depth is 298 or 27. However, when the read depth is up to 3000 level, the sensitivity is  surprisingly decreasing. The sensitivity is 0.71 when coverage is 30590.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He_________________     Date: 09/30/13


Witnessed by: ________________________
