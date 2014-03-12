2013-09-23 SNP calling using varscan2 in somatic mode
==============================

Purpose
------------
To compare varscan2(somatic mode) with rvd2 across coverage depths and dilutions.

Conclusions
-----------------
Run `SNP_varscan2_somatic.sh`to get the final vcf files.

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


The overall experiment is implemented in `SNP_varscan2_somatic.sh`.


Experimental Protocol
---------------------------
**Approach 1.**`SNP_varscan2_somatic1.sh`  **Failed.**

1. mpileup 6 control bam files to `control.pileup`
2. mpileup 6 case bam files to `case100_0.pileup`
3. Feed `control.pileup` and `case100_0.pileup` to varscan somatic

**Result:**

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
	397 positions in tumor
	397 positions shared in normal
	0 had sufficient coverage for comparison
	0 were called Reference
	0 were mixed SNP-indel calls and filtered
	0 were called Germline
	0 were called LOH
	0 were called Somatic
	0 were called Unknown
	0 were called Variant



**Approach 2.**`SNP_varscan2_somatic2.sh`  **Failed.**

1. mpileup 6 control bam files and 6 100.0% bam files to one `control_case100_0.pileup`

1. Feed the single pileup file to varscan somatic. 
 
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
	Reading input from pileup2/10000/control_case100_0.pileup
	Reading mpileup input...
	400 positions in mpileup file
	372 had sufficient coverage for comparison
	372 were called Reference
	0 were mixed SNP-indel calls and filtered
	0 were called Germline
	0 were called LOH
	0 were called Somatic
	0 were called Unknown
	0 were called Variant


The result seems reasonable, with 21 called somatic. **However**, if I ran this program on other dilutions bam files, the result returned are exactly the same, which is definitely what we expected.

**Approach 3.**`SNP_varscan2_somatic3.sh` **Succeed**

1. Use samtools to **merge** 6 bam files(depth 10000, dilution 100.0%) for same dilution to `case100_0.bam` first, mpileup `case100_0.bam` to `case100_0.pileup`

1. Do the same processing to 6 control bam files and get `control.bam` and `control.pileup`

1. Feed `control.pileup` and `case100_0.pileup` to varscan somatic.

**Result:**

	Normal Pileup: pileup3/10000/control.pileup
	Tumor Pileup: pileup3/10000/case100_0.pileup
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
	397 positions in tumor
	397 positions shared in normal
	392 had sufficient coverage for comparison
	375 were called Reference
	0 were mixed SNP-indel calls and filtered
	0 were called Germline
	1 were called LOH
	16 were called Somatic
	0 were called Unknown
	0 were called Variant




**Approach 4.**`SNP_varscan2_somatic4.sh` **Succeed**

1. Use samtools to merge 6 bam files(depth 10000, dilution 100.0%) for same dilution `case100_0.bam`

1. Do the same processing to 6 control bam files and get `control.bam`

1. mpileup `control.bam` and `case100_0.bam` to `control-case100_0.pileup`

1. Feed `control_case100_0.pileup` to varscan somatic.

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
	Parsing Exception on line:
	1_1_400 387     A       75      .$.$,$,$....,,,,,,,,,,..,..,,...,,,,,,,....,,,,.,,,,,nt,...,,,,,,,,,,,,,,,,,,,, @@@@T]^^^^^^^^^^^^caeeaeefeGdfff^ffdfdfc^df^eeffeBJSe^fcffffffffffLBdcffe`f       0
	7

But from the SNP file generated it succeeded in calling somatics. there should be something wrong in mpileup control and case in to one single .pileup file.


**Approach 5.**`SNP_varscan2_somatic5.sh` **Succeed**

mpile only one sorted control bam file and one sorted case bam file to `control.pileup` and `case100_0.pileup` and then feed the two pileup files into varscan2 somatic.

Data:

	control: ../2013-08-06_Downsample_Read_Depth/bam/10000/20100916_c1_p1.02_ACT.bam
	case: 100% dilution: ../2013-08-06_Downsample_Read_Depth/bam/10000/20100916_c3_p1.07_CGT.bam

Result:

	Normal Pileup: pileup5/10000/control.pileup
	Tumor Pileup: pileup5/10000/case100_0.pileup
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
	396 positions in tumor
	396 positions shared in normal
	371 had sufficient coverage for comparison
	357 were called Reference
	0 were mixed SNP-indel calls and filtered
	0 were called Germline
	0 were called LOH
	14 were called Somatic
	0 were called Unknown
	0 were called Variant


**Approach 6.**`SNP_varscan2_somatic6.sh` **Succeed**

mpile only one sorted control bam file and one sorted case bam file together to `control_case100_0.pileup` the pileup file into varscan2 somatic.

Data:

	control: ../2013-08-06_Downsample_Read_Depth/bam/10000/20100916_c1_p1.02_ACT.bam
	case: 100% dilution: ../2013-08-06_Downsample_Read_Depth/bam/10000/20100916_c3_p1.07_CGT.bam

Result:
	Not showing in correct format.
	Did called VarScan somatic correctly in .snp file.

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
	Reading input from pileup6/10000/control_case100_0.pileup
	Reading mpileup input...
	Parsing Exception on line:
	1_1_400 387     A       5       ..,,,   fdcff   0
	7



**Approach 7.**`SNP_varscan2_somatic7.sh`      **Failed**

Same process as **Approach 5.** but with undownsampled data

Data:

	control: ../../data/Synthetic_BAM_files/20100916_c1_p1.02_ACT.bam
	case 100% dilution: ../../data/Synthetic_BAM_files/20100916_c3_p1.07_CGT.bam

Result:

	Normal Pileup: pileup7/10000/control.pileup
	Tumor Pileup: pileup7/10000/case100_0.pileup
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
	395 had sufficient coverage for comparison
	3 were called Reference
	11 were mixed SNP-indel calls and filtered
	131 were called Germline
	73 were called LOH
	2 were called Somatic
	175 were called Unknown
	0 were called Variant



Output of varscan2 somatic 
--------------------------
Normal/Tumor Input Parsing
In somatic mode, VarScan reads the pileup files from normal and tumor simultaneously. Only positions that are present in both files, and meet the minimum coverage in both files, will be compared. VarScan expects that positions on the same chromosome occur in ascending order, so input should be position-sorted! Please note, this kind of simultaneous parsing gets tricky when there are numerous reference sequences (e.g. unplaced contigs) to which reads from only one sample aligned. VarScan tries to obtain the maximum number of comparisons, even if it means closing and reopening the normal file to try to match contig and position. This can lead to looping errors; please report them if you see a number of warmings about "Resetting normal file..." and (if possible) send us sample pileups. 

Variant Calling and Comparison
At every position where both normal and tumor have sufficient coverage, a comparison is made. First, normal and tumor are called independently using the germline consensus calling functionality. Then, their genotypes are compared by the following algorithm: 

    If tumor matches normal: 
        If tumor and normal match the reference 
            ==> Call Reference 
        Else tumor and normal do not match the reference 
            ==> Call Germline 

    Else tumor does not match normal: 
    Calculate significance of allele frequency difference by Fisher's Exact Test 
        If difference is significant (p-value < threshold):
            If normal matches reference
                ==> Call Somatic
            Else If normal is heterozygous
                ==> Call LOH
            Else normal and tumor are variant, but different
                ==> Call IndelFilter or Unknown
        Else difference is not significant:
        Combined tumor and normal read counts for each allele. Recalculate p-value.
            ==> Call Germline 

So what we are expecting should be 'Somatic', 'LOH' and 'IndelFilter' or 'Unknowncalled' in variant positions.
Results
-----------
Approaches which didn't merge bam files:

**Failed.** **Approach 1** is not correct since it returns the message **'0 had sufficient coverage for comparison'**, which means for some reason the pileup files are not correctly read by VarScan somatic. 

**Failed.** **Approach 2** called no variants. All positions are called as reference. 
****************************

Approaches which did merge 6 bam files (replicates):

**Succeed** **Approach 3** 1 were called LOH, 16 were called Somatic

**Succeed** **Approach 4** 1 were called LOH, 16 were called Somatic

Same results. From the .snp files actually 1 were called LOH,14 were called Sometic
****************************

Approaches where one single bam file for control and one single bam file for case

**Succeed** **Approach 5** 14 were called Somatic

**Succeed** **Approach 6** 1 14 were called Somatic

Most accurate.
****************************
undownsampled dataset:

**Failed** **Approach 7** 
****************************

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________     Date: _____________________


Witnessed by: ________________________
