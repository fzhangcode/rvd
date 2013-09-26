2013-09-23 SNP calling using varscan2 in somatic mode
==============================

Purpose
------------
To compare varscan2(somatic mode) with rvd2 across coverage depths and dilutions.

Conclusions
-----------------

Background
----------------

Materials and Equipment
------------------------------
rvd2/bin/varscan.v2.3.4.jar somatic

Series of experiment on different approaches were done as : `SNP_varscan2_somatic2.sh`, `SNP_varscan2_somatic3.sh`, `SNP_varscan2_somatic4.sh`. The data used in these experiments are:

	DOWNDIR=../2013-08-06_Downsample_Read_Depth/bam/10000
	
	control data: 
	SORTCONTROL=$(ls $DOWNDIR/20100916_c1_p?.02*.sorted.bam $DOWNDIR/20100916_c1_p?.04*.sorted.bam $DOWNDIR/20100916_c1_p?.05*.sorted.bam)
	case data: 
	SORTBAM100_0=$(ls $DOWNDIR/20100916_c3_p?.07*.sorted.bam  $DOWNDIR/20100916_c3_p?.12*.sorted.bam $DOWNDIR/20100916_c3_p?.14*sorted.bam)


The overall experiment is implemented in `SNP_varscan2_somatic.sh`, which is not yet working right now.


Experimental Protocol
---------------------------
**Approach 1.**`SNP_varscan2_somatic1.sh` 

1. mpileup 6 control bam files to `control.pileup`
2. mpileup 6 case bam files to `case100_0.pileup`
3. Pass `control.pileup` and `case100_0.pileup` to varscan somatic

**Result:**

	Pass the single pileup file to varscan somatic. 
	Normal Pileup: pileup1/10000/control.pileup
	Tumor Pileup: pileup1/10000/case100_0.pileup
	Min coverage:   8x for Normal, 6x for Tumor
	Min reads2:     2
	Min strands2:   1
	Min var freq:   1.0E-5
	Min freq for hom:       0.75
	Normal purity:  1.0
	Tumor purity:   1.0
	Min avg qual:   15
	P-value thresh: 0.99
	Somatic p-value:        0.05

	400 positions in tumor
	400 positions shared in normal
	**0 had sufficient coverage for comparison**
	0 were called Reference
	0 were mixed SNP-indel calls and filtered
	0 were called Germline
	0 were called LOH
	0 were called Somatic
	0 were called Unknown
	0 were called Variant


**Approach 2.**`SNP_varscan2_somatic2.sh` 

1. mpileup 6 control bam files and 6 100.0% bam files to one `control_case100_0.pileup`

1. Pass the single pileup file to varscan somatic. 
 
**Result:**

	Min coverage:   8x for Normal, 6x for Tumor
	Min reads2:     2
	Min strands2:   1
	
	400 positions in mpileup file
	400 had sufficient coverage for comparison
	4 were called Reference
	0 were mixed SNP-indel calls and filtered
	261 were called Germline
	0 were called LOH
	21 were called Somatic
	114 were called Unknown
	0 were called Variant

The result seems reasonable, with 21 called somatic. **However**, if I ran this program on other dilutions bam files, the result returned are exactly the same, which is definitely what we expected.

**Approach 3.**`SNP_varscan2_somatic3.sh`

1. Use samtools to **merge** 6 bam files(depth 10000, dilution 100.0%) for same dilution to `case100_0.bam` first, mpileup `case100_0.bam` to `case100_0.pileup`

1. Do the same processing to 6 control bam files and get `control.bam` and `control.pileup`

1. Pass `control.pileup` and `case100_0.pileup` to varscan somatic.

**Result:**

	Normal Pileup: pileup3/10000/control.pileup
	Tumor Pileup: pileup3/10000/case100_0.pileup
	Min coverage:   8x for Normal, 6x for Tumor
	Min reads2:     2

	400 positions in tumor
	400 positions shared in normal
	400 had sufficient coverage for comparison
	6 were called Reference
	8 were mixed SNP-indel calls and filtered
	176 were called Germline
	36 were called LOH
	0 were called Somatic
	174 were called Unknown
	0 were called Variant



**Approach 4.**`SNP_varscan2_somatic4.sh`

1. Use samtools to merge 6 bam files(depth 10000, dilution 100.0%) for same dilution `case100_0.bam`

1. Do the same processing to 6 control bam files and get `control.bam`

1. mpileup `control.bam` and `case100_0.bam` to `control-case100_0.pileup`

1. Pass `control_case100_0.pileup` to varscan somatic.

**Result:**

	Min coverage:   8x for Normal, 6x for Tumor
	Min reads2:     2
	Min strands2:   1
	Min var freq:   1.0E-5
	Min freq for hom:       0.75
	Normal purity:  1.0
	Tumor purity:   1.0
	Min avg qual:   15
	P-value thresh: 0.99
	Somatic p-value:        0.05
	Reading input from pileup4/10000/control_case100_0.pileup
	Reading mpileup input...

	400 positions in mpileup file
	400 had sufficient coverage for comparison
	6 were called Reference
	8 were mixed SNP-indel calls and filtered
	176 were called Germline
	36 were called LOH
	0 were called Somatic
	174 were called Unknown
	0 were called Variant



Results
-----------
**Approach 1** is not correct since it returns the message **'0 had sufficient coverage for comparison'**, which means for some reason the pileup files are not correctly read by VarScan somatic.

**Approach 2** called 21 somatic. which seems reasonable at the first sight. **However**, if we ran this program on other dilutions bam files, the result returned are exactly the same, which is definitely what we expected. I doubt that the case pileup are not read by VarScan somatic so there is no difference between dataset with different dilution rate.


**Approach 3** and **Approach 4** got the same result. However, no positions are called as variant or somatic. 

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________     Date: _____________________


Witnessed by: ________________________
