2013-10-01 SNP calling using strelka
==============================

Purpose
------------
To compare strelka with rvd2 on variant calling

Conclusions
-----------------
This approach is not able to run on our dataset for certain reason.

Background
----------------
Strelka is an analysis package designed to detect somatic SNVs and small indels from the aligned sequencing reads of matched tumor-normal samples. Webpage
[Strelka somatic variant caller - Google Sites](https://sites.google.com/site/strelkasomaticvariantcaller/) contains all the information about Strelka.


Materials and Equipment
------------------------------



Experimental Protocol
---------------------------
A sample experiment is implemented in `SNP_strelka.sh`

Results
-----------

I got the essentially same bug as when I ran Virmid.

So the error returned was:  

	ERROR: BAM headers and reference fasta disagree on chromosome: '2_1_10'
	'BAM headers and reference fasta disagree on chromosome: \'2_1...'

I used samtools to view the headers of the sorted bam file and got:

	@HD     VN:1.4  SO:coordinate
	@SQ     SN:1_1_400      LN:400
	@SQ     SN:2_1_10       LN:10
	@PG     ID:bwa  PN:bwa  VN:0.5.9-r16

Why there is a chromosome 2_1_10 in the bam file header?


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____Yuting He____________     Date: _______09/30/13______________


Witnessed by: ________________________
