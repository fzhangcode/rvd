2014-01-30 Variants Calling on Yeast Data Experiment
=============

Purpose: 
---------------------
Call mutations by testing generation 007 vs generation 133 from E3 on GSY1135_Chr4. 

Call mutations by testing generation 007 vs generation 300 from E3 on GSY1135_Chr4.

Try to find the gene MTH1&IRA1, MNN4, HXT617 by testing 1000 bases (Chr4:1014401:1015702).

Materials:
--------
	control: generation 007 (SRR515487)
	case1: generation 133 (SRR519064) 
	case2: generation 322 (SRR519085)
~~

	add the executables to your path:
	Echo $PATH
	Export PATH=$PATH:/home/fzhang/Research/rvd2/results/2013-12-20_Yeast_data_processing/bwa-0.7.5a)
~~
	
	To be more clear, below is the content for the folder:
	GSY1135_Chr4: the reference sequence for GSY1135 on Chr4
	S288C_Reference: all the reference of S288C from Chr1 to Chr16
	samples_for_Chr4:samples of experiment for Chr4, including the data file of FASTQ, SAM, BAM, Pileup, DC and Intermediate folder. 

Conclusions
-----------------
1. RVD3 with lognormal prior didn't call any mutations.

2. The mutation of MTH1 doesn't show up in the depth chart file, so there must be some step wrong when processing the data from SRA format to the depth chart files.

Background
-----------------
From Kvitek's paper, they called mutations by three experiments. In our research, we want to find MTH1 which is in Chr4 in Table_S2 by kvitek, 2013:

position - Ref -	Alt

1014583 - G	-  A

1014691  - G  - A

1014981 - A - T

1015077 - G - T

So we choose the positions between Chr4:1014401:1015702 (1000 postitions) on Chromosome 4 for testing to find mutations.

1.	Reference sequence for GSY1135
-------
Download GSY1135 from [27] 

	http://sra.dnanexus.com/studies/SRP002895/runs (SRA020606.1 -> SRR063399)
	Save it to Y:/freeze/kvitek2013/GSY1135_Chr4

Map GSY1135 to S288C (Chr4) using BWA (GSY1135: reads; S288C: references)

	Indexd fastq: bwa index S288C_Reference/Chr4_Reference/Chr04.fa

	Align indexd fastq:    bwa aln –B 6 –t 8 –q 10 S288C_Reference/Chr4_Reference/Chr4.fa GSY1135_Chr4/SRR063399.fastq > GSY1135_Chr4/SRR063399.sai

	bwa samse S288C_Reference/Chr4_Reference/Chr4.fa GSY1135_Chr4/SRR063399.sai GSY1135_Chr4/SRR063399.fastq > GSY1135_Chr4/SRR063399.sam

	Convert SAM to BAM:	samtools view –bht Chr4_Reference/Chr4.fa –o ./GSY1135_Chr4/SRR063399.bam ./GSY1135_Chr4/SRR063399.sam

	Sort the BAM file:	 samtools sort ./SRR063399.bam ./SRR063399.sort

	Index the sorted BAM file:	 samtools index ./SRR063399.sort.bam		


Call SNPs using GATK UnifiedGenotyper 
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

	Download picard-tools-1.104.zip from:	 http://sourceforge.net/projects/picard/files/

	Create Dictionary:	java -jar CreateSequenceDictionary.jar R=/home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr4_Reference/Chr04.fa o=/home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr4_Reference/Chr04.dict

	Add head:	Java –jar AddOrReplaceReadGroups.jar I=SRR063399.sort.bam O=./headSRR063399.sort.bam LB=whatever PL=illumine PU=whatever SM=whatever

	Index:	Samtools index ./headSRR063399.sort.bam

	Call SNPs:		Java –jar GenomeAnalysisTK.jar –R /home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr4_Reference/Chr4.fa –T UnifiedGenotyper –I /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr4/headSRR063399.sort.bam –o SRR063399.vcf –stand_call_conf 50.0 –stand_emit_conf 10.0
	

Create a FASTA GSY1135 reference using GATK FastaAlternative (GSY1135\_Chr10)
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html 


	java –Xmx2g –jar GenomeAnalysisTK.jar –R /home/pjflaherty/flahertylab/freeze/kvitek2013/S288C_Reference/Chr4_Reference/Chr4.fa –T FastaAlternateReferenceMaker  –O  /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr4/GSY1135_Chr4.fasta –-variant /home/pjflaherty/flahertylab/freeze/kvitek2013/GSY1135_Chr4/SRR063399.vcf

	Create Dictionary:  Java –jar CreateSequeneDictionary.jar R=GSY1135_Chr4.fasta O=GSY1135_Chr4.dict

2.	Download SRAs from Dan’s Paper
-------	
	http://sra.dnanexus.com/studies/SRP013879/runs
	Save the SRA files to the samples_for_Chr4

3.	Convert SRA0 and SRA100 to FASTQ format
-------
	Download the SRA toolkit from http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software 
	Decompress: Tar –zxvf sratoolkit.2.3.4-2-centos_linux64.tar.gz
	Convert SRA to fastq: ./sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump /home/pjflaherty/flahertylab/freeze/kvitek2013/samples_for_Chr4/SRA/*.sra


4.	Remove WT population using FASTX Barcode Splitter
--------
	http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage

	Under the directory: samples_for_Chr10/Intermediate
	cat /home/pjfhaherty/flahertylab/freeze/kvitek2013/samples_for_Chr4/FASTQ/SRR515487.fastq | bin/fastx_barcode_splitter.pl --bcfile /home/pjflaherty/flahertylab/freeze/kvitek2013/samples_for_Chr4/FASTQ/MyBarCodes.txt --bol --mismatches 2 --prefix ./bla_ --suffix ".txt"

[Past. Trim Nextera tag using Cutadapt]	
------	
	https://pypi.python.org/pypi/cutadapt
	cutadapt –a CAAGCAGAAGACGGCATACGAGATNNNNNNCGGTCTGCCTTGCCAGCCCGCTCAG –m 15 Gen007_unmatched.fastq > Gen007_trimed.fastq 

5. Cut pair ends: 
---
	http://www.biostars.org/p/19446/
	awk 'NR%2==1 { print $0 "/1" } ; NR%2==0 { print substr($0,0,length($0)/2) }' /samples_for_Chr4/Intermediate/gen007_unmatched.txt > /samples_for_Chr4/Intermediate/gen007_test1.fastq
	awk 'NR%2==1 { print $0 "/2" } ; NR%2==0 { print substr($0,length($0)/2) }' /samples_for_Chr4/Intermediate/gen007_unmatched.fastq > /samples_for_Chr4/Intermediate/gen007_test2.fastq

6. Map FASTQ0 and FASTQ100 to GSY1135 reference using BWA.
----

    Index reference: bwa index ./GSY1135_Chr4/GSY1135_Chr4.fasta
    samtool faidx ./GSY1135_Chr4/GSY1135_Chr4.fasta
    
	[bwa aln -t 8 -I -q 10 ./GSY1135_Chr10/GSY1135_Chr10.fasta Intermediate/gen007_test1.fastq > Intermediate/gen007_test1.sai
	bwa aln -t 8 -I -q 10 ./GSY1135_Chr10/GSY1135_Chr10.fasta Intermediate/gen007_test2.fastq > Intermediate/gen007_test2.sai
	]
	
	bwa mem ./GSY1135_Chr4/GSY1135_Chr4.fasta ./samples_for_sample4/Intermediate/gen007_test1.fastq ./samples_for_sample4/Intermediate/gen007_test2.fastq > SAM/gen007.sam
	[errors when using sampe]

	Convert SAM to BAM:
	samtools view -b -S -o ./BAM/gen007.bam ./SAM/gen007.sam


7. Make pileup files for FASTA0 and FASTQ100 using Samtools
---
	Add head: java -jar AddOrReplaceReadGroups.jar I=../../../../../../../freeze/kvitek2013/samples_for_Chr4/BAM/gen007.bam o=../../../../../../../freeze/kvitek2013/samples_for_Chr4/BAM/headgen007.bam LB=whatever PL=illumine PU=whatever SM=whatever

	samtools sort BAM/headgen007.bam BAM/headgen007.sort

	samtools index /headgen007.sort.bam

    samtools mpileup -C 50 -d 100000 -f ./GSY1135_Chr4/GSY1135_Chr4.fasta samples_for_Chr4/BAM/headgen0074.sort.bam > samples_for_Chr4/Pileup/gen007.pileup


8.	Make depthchart files using pileup2dc
------
	cd src/pileup2dc
	gcc -o pileup2dc main.c
	mv pileup2dc ../../bin/
	PATH=$PATH:/home/fzhang/Research/rvd2/bin
	cd fzhang/Research/rvd2/results/2014-01-05-pileup/
	python test_pileup2dc.py      
	[error: -> Segmentation fault (core dumped) or the output file is empty!]

    Problem solved: Copy main.c to 2014-01-05-pileup, and revise the "allocate space for the depth chart" with pile_t pile[1], so the segmentation fault won't be happened.
	cd results/2014-01-05-pileup/backup/ 
    gcc -o pileup2dc main.c
    ./pileup2dc ../../../../../../freeze/kvitek2013/samples_for_Chr4/Pileup/gen007.mpileup gen007.dc      

9.	Estimate models for case + control
-----
	Under folder results/2014-01-30_Yeast_Chr4_Jeffreys_test
	Get hdf5 file: python runall.py
	Hypotheses test: python test_rvd29.py
	Summary the results: python character.py

10.	Test for variants between case + control -> vcf files
-----


###Note1: 
	http://biobits.org/samtools_primer.html

###Note2: linux
	control + A + D
    screen -S/ -ls/ -r
    ls -lh 

###Note3: spiked fresh of GSY1135
	#matched reads/(#matched+#unmatched)=spiked fresh of GSY1135


Prepared by: _________ _Fan Zhang______ Date: ____________12/20/2013____________


Witnessed by: ________________________