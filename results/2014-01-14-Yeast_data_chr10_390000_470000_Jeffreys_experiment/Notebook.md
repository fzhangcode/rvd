2014-02-04 Variants Calling on Yeast Data Experiment 
==============================

Purpose
------------
Call mutations on the positions from ChrX:425157 to 431237 (gene CYR1) by testing generation 007 vs generation 133 from experiment 1 (Kvitek,2013).


Conclusion
--------
Since we call different mutations compared with Kvitek's work, there may be some step wrong in the preprocessing from SRA format to the depth chart file, or the reference generation process.


Background
-----------------
From Kvitek's paper, they called mutations by three experiments. In our research, we want to call mutations from ChrX:425157 to 431237 by RVD3 model with variations in priors function. So we choose the positions between ChrX:390000 and 470000 (8000 positions) on Chromosome 10 for experiment.


Materials and Equipment
------------------------------
control data: SRR515969 (generation007)

case data: SRR519089 (generation133)


Experimental Protocol
---------------------------
1.	Generate reference sequence for GSY1135
-------
Download GSY1135 

	http://sra.dnanexus.com/studies/SRP002895/runs (SRA020606.1 -> SRR063399)
	Save it to Y:/freeze/kvitek2013/GSY1135_Chr10/

Map GSY1135 to S288C (Chr10) using BWA (GSY1135: reads; S288C: references)

	Indexd fastq: "bwa index S288C_Reference/Chr10_Reference/Chr10.fa"

	Align indexd fastq: "bwa aln –B 6 –t 8 –q 10 S288C_Reference/Chr10_Reference/Chr10.fa GSY1135_Chr10/SRR063399.fastq > GSY1135_Chr10/SRR063399.sai"

	Generate alignments in the SAM format given signle-end reads: 
	"bwa samse S288C_Reference/Chr10_Reference/Chr10.fa GSY1135_Chr10/SRR063399.sai GSY1135_Chr10/SRR063399.fastq > GSY1135_Chr10/SRR063399.sam"

	Convert SAM to BAM:	"samtools view –bht Chr10_Reference/Chr10.fa –o ./GSY1135_Chr10/SRR063399.bam ./GSY1135_Chr10/SRR063399.sam"

	Sort the BAM file:	 "samtools sort ./SRR063399.bam ./SRR063399.sort"

	Index the sorted BAM file: "samtools index ./SRR063399.sort.bam"		


Call SNPs using GATK UnifiedGenotyper 
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

	Download picard-tools-1.104.zip from:	 http://sourceforge.net/projects/picard/files/

	Create Dictionary:	java -jar CreateSequenceDictionary.jar R=/home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr10_Reference/Chr10.fa o=/home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr10_Reference/Chr10.dict

	Add head:	Java –jar AddOrReplaceReadGroups.jar I=SRR063399.sort.bam O=./headSRR063399.sort.bam LB=whatever PL=illumine PU=whatever SM=whatever

	Index:	Samtools index ./headSRR063399.sort.bam

	Call SNPs:	Java –jar GenomeAnalysisTK.jar –R /home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr10_Reference/Chr10.fa –T UnifiedGenotyper –I /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr10/headSRR063399.sort.bam –o SRR063399.vcf –stand_call_conf 50.0 –stand_emit_conf 10.0
	

Create a FASTA GSY1135 reference using GATK FastaAlternative (GSY1135\_Chr10)
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html 


	Generate reference sequence: "java –Xmx2g –jar GenomeAnalysisTK.jar –R /home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr10_Reference/Chr10.fa –T FastaAlternateReferenceMaker  –O  /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr10/GSY1135_Chr10.fasta –-variant /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr10/SRR063399.vcf"

	Create Dictionary:  "Java –jar CreateSequeneDictionary.jar R=GSY1135_Chr10.fasta O=GSY1135_Chr10.dict"

2.	Download SRAs from Kvitek’s Paper
-------	
	http://sra.dnanexus.com/studies/SRP013879/runs
	Save SRR515969.sra and SRR519089.sra files to the samples_for_Chr10 folder

3.	Convert SRA to FASTQ format
-------
	Download the SRA toolkit from http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software 
	Decompress: "Tar –zxvf sratoolkit.2.3.4-2-centos_linux64.tar.gz"
	Convert SRA to fastq: "./sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump /home/pjflaherty/flahertylab/freeze/kvitek2013/samples_for_Chr10/SRA/*.sra"


4.	Remove WT population using FASTX Barcode Splitter
--------
	http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage
	Sample was extracted with the exact tag ATCTCG:
	"cat /home/pjfhaherty/flahertylab/freeze/kvitek2013/samples_for_Chr10/FASTQ/SRR515969.fastq | bin/fastx_barcode_splitter.pl --bcfile /home/pjflaherty/flahertylab/freeze/kvitek2013/samples_for_Chr10/FASTQ/MyBarCodes.txt --bol --mismatches 2 --prefix ./bla_ --suffix ".txt""


5. Cut pair-end reads
---
	http://www.biostars.org/p/19446/ 
	Both pair-end reads are joined in one FASTQ. I need to split the file in two separated FASTQ pair-end files. Sequence identifiers will have /1 or /2 appended for the split left-hand and right-hand reads, respectively.
	For the left read: "awk 'NR%2==1 { print $0 "/1" } ; NR%2==0 { print substr($0,0,length($0)/2) }' /samples_for_Chr10/Intermediate/gen007_unmatched.fastq > /samples_for_Chr10/Intermediate/gen007_left.fastq"
	For the right read: "awk 'NR%2==1 { print $0 "/2" } ; NR%2==0 { print substr($0,length($0)/2) }' /samples_for_Chr10/Intermediate/gen007_unmatched.fastq > /samples_for_Chr10/Intermediate/gen007_right.fastq"

6. Map FASTQ files to GSY1135 reference using BWA and get BAM files.
----

    Index reference: "bwa index ./GSY1135_Chr10/GSY1135_Chr10.fasta"
    "samtool faidx ./GSY1135_Chr10/GSY1135_Chr10.fasta"
    
	Alignment: "bwa mem ./GSY1135_Chr10/GSY1135_Chr10.fasta ./samples_for_Chr10/Intermediate/gen007_left.fastq ./samples_for_Chr10/Intermediate/gen007_right.fastq > SAM/gen007.sam"

	Convert SAM to BAM:	"samtools view -b -S -o ./BAM/gen007.bam ./SAM/gen007.sam"
	
	Add head to the BAM file: "java -jar AddOrReplaceReadGroups.jar I=../../../../../../../freeze/kvitek2013/samples_for_Chr10/BAM/gen007.bam o=../../../../../../../freeze/kvitek2013/samples_for_Chr10/BAM/headgen007.bam LB=whatever PL=illumine PU=whatever SM=whatever"

	Sort the BAM file: "samtools sort BAM/headgen007.bam BAM/headgen007.sort"

	Index: "samtools index /headgen007.sort.bam"

7. Make pileup files using Samtools
---


    Pileup: "samtools mpileup -C 50 -d 100000 -f ./GSY1135_Chr10/GSY1135_Chr10.fasta samples_for_Chr10/BAM/headgen007.sort.bam > samples_for_Chr10/Pileup/gen007.pileup"


8.	Make depthchart files using pileup2dc
------
	"cd results/2014-01-05-pileup/backup/" 
    "gcc -o pileup2dc main.c"
    "./pileup2dc ../../../../../../freeze/kvitek2013/samples_for_Chr10/Pileup/gen007.mpileup gen007.dc"      

9.	Estimate models for case + control
-----
	"cd results/2014-01-30_Yeast_Chr10_Jeffreys_test"
	Get hdf5 file: "python runall.py"

10.	Test for variants between case + control -> vcf files
---
	Hypotheses test: "python test_rvd29.py"
	Summary the results: "python character.py"

Results
-----------
1. Called mutations:
 
	From the generated VCF files, we call 27 variants from ChrX:390000 to 470000 by RVD3 with lognormal prior. And all these variants are covered in the called 59 variants by RVD3 with Jeffreys prior except one (ChrX:453972).

2. One mutation in gene CYR1 range: 

	Among the mutations we called, one mutant position (ChrX:428201, Ref:A->Alt:C) is located in the CYR1 range (ChrX:425157 to 431237).
 
3. Compare with Kvitek's results:
	
	We failed to call mutation on the position ChrX:428286, which is called as a mutation (Ref:A->Alt:C) by Kvitek (see Table_S2). This mutation doesn't show up in our depth chart file either. 


Archived Samples
-------------------------
All the original data and archived samples are stored in: freeze/kvitek2013/

Archived Computer Data
------------------------------


Prepared by: _______Fan Zhang____________ Date: ______2014.02.04_______________


Witnessed by: ________________________
