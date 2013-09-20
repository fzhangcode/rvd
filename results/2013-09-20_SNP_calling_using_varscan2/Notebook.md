2013-09-20 SNP calling using varscan2
==============================

Purpose
------------
To compare varscan2 with rvd2 across coverage depths and dilutions.

Conclusions
-----------------

Background
----------------

Materials and Equipment
------------------------------
rvd2/bin/varscan.v2.3.4.jar


Experimental Protocol
---------------------------

To process 10x downsampled data from synthetic DNA and a 0.1% dilution. Run this.

samtools mpileup -d 100000 -f ../../data/Synthetic_BAM_files/plasmid.fa ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.07_CGT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.12_GTT.sorted.bam ../2013-08-06_Downsample_Read_Depth/bam/10/20100916_c1_p1.14_TCT.sorted.bam > case0_1.pileup

java -jar ../../bin/VarScan.v2.3.4.jar mpileup2snp case0_1.pileup --output-vcf 1 > case0_1.vcf

Results
-----------


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________     Date: _____________________


Witnessed by: ________________________
