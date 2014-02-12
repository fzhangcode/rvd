2014-2-10 Variants Calling on HCC1187 Data 
=============

Purpose: 
---------------------
Call mutations by RVD3 on HCC1187 data.

Conclusions
-----------------
The results are consistent with Yuting's results.


Background
-----------------
The HCC1187 dataset is a well-recognized baseline dataset from
Illumina for evaluating sequence analysis algorithms (Newman
et al., 2013; Howarth et al., 2011, 2007). The HCC1187 cell line
was derived from epithelial cells from primary breast tissue from a
41 y/o adult with TNM stage IIA primary ductal carcinoma. The
estimated tumor purity was reported to be 0.8. Matched normal
cells were derived from lymphoblastoid cells from peripheral blood.
Sequencing libraries were prepared according to the protocol described
in the original technical report (Allen, 2013). The raw FASTQ
read files were aligned to hg19 using the Isaac aligner to generate
BAM files (Raczy et al., 2013). The aligned data had an average
read depth of 40x for the normal sample and 90x for the tumor
sample with about 96% coverage with 10 or more reads. We used
samtools mpileup to generate pileup files using hg19 as reference
sequence (Navin et al., 2010).

Materials and Equipment
------------------------------
Run the depth chart data directly from: 

	/yhe2/Research/rvd2/data/HCC1187_PAXIP1_genome_depth_chart_600/. 

	control: HCC1187BL_S1.dc
	case: HCC1187C_S1.dc 

The generated hdf5 and vcf files are saved in the current directory: 

	/fzhang/Research/rvd2/results/2014-02-07_RVD3_test_on_HCC1187/.


Experimental Protocol
---------------------------
Generate Hdf5:
	
	 Run runall\_case.py and runall\_control.py simultaneously 

Generate Vcf: 

	Run test_rvd29.py  

Results
-----------------
The RVD3 with Jeffreys prior called 3 mutations, which are the same with the results from RVD3 with lognormal prior:

CHROM	POS		ID	REF	ALT	

7	154754371	.	T	C	

7	154758813	.	G	A	

7	154760439	.	A	C	

chr7:154754371T>C, chr7:154758813G>A, chr7:154760439A>C
are significantly different in tumor and normal sample.

Positions chr7:154754371T>C and
chr7:154758813G>A are called loss-of-heterozygosity events
	
Archived Samples
-------------------------


Archived Computer Data
------------------------------



Prepared by: _________ _Fan Zhang______ Date: ____________2/10/2014____________


Witnessed by: ________________________