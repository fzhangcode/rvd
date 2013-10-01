2013-09-20 SNP calling using varscan2
==============================

Purpose
------------
To compare varscan2 with rvd2 across coverage depths and dilutions.

Conclusions
-----------------
VarScan2 mpileup2snp is only able to call 100% snps.

Background
----------------
Please go to [VarScan2 mpileup2snp manual page](http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp) for more information.

Materials and Equipment
------------------------------
rvd2/bin/varscan.v2.3.4.jar


Experimental Protocol
---------------------------

To process 10x downsampled data from synthetic DNA and a 0.1% dilution. Run this.

This generates case pileups from pair 1 data
`samtools mpileup -d 100000 -C 50 -f ../../data/Synthetic_BAM_files/plasmid.fa ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.07_CGT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.12_GTT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.14_TCT.sorted.bam > case-10x-0_1.pileup`

This generates combine control pileups from pair1 data
`samtools mpileup -d 100000 -C 50 -f ../../data/Synthetic_BAM_files/plasmid.fa ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.02_ACT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.04_ATT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.05_CAT.sorted.bam > control-10x.pileup`

This uses the somatic version of varscan
`java -jar ../../bin/VarScan.v2.3.4.jar somatic control-10x.pileup case-10x-0_1.pileup blah --tumor-purity 0.001`

This only looks for 100% snps
`java -jar ../../bin/VarScan.v2.3.4.jar mpileup2snp case0_1.pileup --output-vcf 1 > case0_1.vcf`


The `SNP_varscan2` applies `mpileup2snp` across all datasets for variant calling.

Results
-----------
The program generates vcf files which saves variant calling results. They are available in the same directory.

Using the characer.py in `../2013-09-19_operating_characteristics` it is possible to generate a summarizing table showing the sensitivity/specificity for variant calling using VarScan2 mpileup2snp. A screenshot is shown as follows:

![](http://i.imgur.com/KMgF7rx.png)

Figure 1. Sensitivity/specificity achieved using VarScan2 mpileup2snp for variant calling.


It can be seen that VarScan2 mpileup2snp is only able to call 100% snps. The sensitity and specificity achieve is 1.00 and 1.00 across all coverage options. 



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Yuting He________     Date: _________09/30/13____________


Witnessed by: ________________________
