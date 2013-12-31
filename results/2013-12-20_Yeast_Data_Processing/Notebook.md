Variants Calling on Yeast Data Experiment Outline - 10 Steps
=============

Goal: What are the mutations?  At generation 0 vs generation 100? Choose GSY1135_Chr10. Look for where is CYR1 in E1)
------
Materials:
--------
	SRA for gen0: SRA0 = control (E1_gen007=SRR515969)
	SRA for gen100: SRA100 = case (E1_gen133=SRR519089)
~~

	add the executables to your path:
	Echo $PATH
	Export PATH=$PATH:/home/fzhang/Research/2013_research_plan/Yeast_experiment/bwa-0.7.5a)

1.	Reference sequence for GSY1135
-------
Download GSY1135 from [27] 

	http://sra.dnanexus.com/studies/SRP002895/runs (SRA020606.1 -> SRR063399)

Map GSY1135 to S288C (Chr10) using BWA (GSY1135: reads; S288C: references)

	Align indexd fastq:    bwa aln –B 6 –t 8 –I –q 10 ./Chr10.fa ./Reference_GSY1135/SRR063399.fastq > ./SAI_files/SRR063399.sai

	bwa samse ./Chr10.fa ./SAI_files/SRR063399.sai ./Reference_GSY1135/SRR063399.fastq > ./Reference_GSY1135/ SRR063399.sam

	Convert SAM to BAM:	samtools view –bt ./Chr10.fa –o ./Reference_GSY1135/SRR063399bam ./ Reference_GSY1135/SRR063399.sam

	Sort the BAM file:	 samtools sort ./SRR063399.bam ./SRR063399.sort

	Index the sorted BAM file:	 samtools index ./SRR063399.sort.bam		rename the reference -> GSY1135_Chr10

Call SNPs using GATK UnifiedGenotyper 
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

	Download picard-tools-1.104.zip from:	 http://sourceforge.net/projects/picard/files/

	Create Dictionary:	Java –jar CreateSequeneDictionary.jar R=Chr10.fa O=Chr10.dict

	Add head:	Java –jar AddOrReplaceGroups.jar I=SRR063399.sort.bam O=./headSRR063399.sort.bam LB=whatever PL=illumine PU=whatever SM=whatever

	Index:	Samtools index ./headSRR063399.sort.bam

	Call SNPs:		Java –jar GenomeAnalysisTK.jar –R Chr100.fa –T UnifiedGenotyper –I  ./headSRR063399.sort.bam –o SRR063399.vcf –stand_call_conf 50.0 –stand_emit_conf 10.0

Create a FASTA GSY1135 reference using GATK FastaAlternative (GSY1135\_Chr10)
http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_fasta_FastaAlternateReferenceMaker.html 

	java –Xmx2g –jar GenomeAnalysisTK.jar –R Chr10.fa –T FastaAlternateReferenceMaker  –O  GSY1135_Chr10.fasta –variant ./Reference_GSY1135/SRR063399.vcf

	Create Dictionary:  Java –jar CreateSequeneDictionary.jar R=GSY1135_Chr10.fasta O=GSY1135_Chr10.dict

2.	Download SRAs from Dan’s Paper
-------	
	http://sra.dnanexus.com/studies/SRP013879/runs

3.	Convert SRA0 and SRA100 to FASTQ format
-------
	Download the SRA toolkit from http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software 
	Decompress: Tar –zxvf sratoolkit.2.3.4-2-centos_linux64.tar.gz
	Convert SRA to fasta: ./sratoolkit.2.3.4-2-centos_linux64/bin/fastq-dump /home/pjflaherty/flahertylab/freeze/baker_yeast/*.sra

4.	Remove WT population using FASTX Barcode Splitter
--------
	http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage

	cat Gen007.fastq| /usr/local/bin/fastx_barcode_splitter.pl --bcfile mybarcodes.txt --bol --mismatches 2 \ --prefix /tmp/bla_ --suffix ".txt

5.	Trim Nextera tag using Cutadapt	
------	
	https://pypi.python.org/pypi/cutadapt
	cutadapt –a CAAGCAGAAGACGGCATACGAGATNNNNNNCGGTCTGCCTTGCCAGCCCGCTCAG –m 15 Gen007_unmatched.fastq > Gen007_trimed.fastq 

6. Map FASTQ0 and FASTQ100 to GSY1135 reference using BWA. 
--------
   (Quality score: 50)	

	Align indexd fastq:	bwa aln –B 6 –t 8 –I –q 10 ./Reference_GSY1135/GSY1135_Chr10.fasta ./BARCODE/Gen007_trimed.fastq > ./SAI_files/Gen007.sai

	bwa samse ./Reference_GSY1135/GSY1135_Chr10.fasta ./SAI_files/Gen007.sai ./BARCODE/Gen007_trimed.fastq > ./Control_Case_files/Gen007.sam

	Convert SAM to BAM:	
	samtools view –bt ./Reference_GSY1135/ GSY1135_Chr10.fasta –o ./Control_Case_files/ Gen007.bam ./Control_Case_files/Gen007.sam
	or 
	samtools view –b –S –o ./Gen007_1.bam ./Gen007.sam

	(Sorted BAM file was created with Picard v1.45 FixMateInformation   http://picard.sourceforge.net/
	java -Xmx4g -jar FixMateInformation.jar INPUT=file.realigned.bam OUTPUT=file.realigned.fixedmateinfo.bam SO=coordinate MAX_RECORDS_IN_RAM=5000000 VALIDATION_STRINGENCY=LENIENT  CREATE_INDEX=true TMP_DIR=$TMPDIR)

7.	Make pileup files for FASTA0 and FASTQ100 using Samtools
------
	http://www.biostars.org/p/63429/#73512

	Add head:   Java –jar AddOrReplaceGroups.jar I=Gen007.bam O=./headGen007.bam LB=whatever PL=illumine PU=whatever SM=whatever

	Sort the BAM file:	 samtools sort -n ./Gen007.bam ./Gen007.sort

	Index:	Samtools index ./Gen007.sort.bam

	samtools mpileup –C 50 -uf ./Reference_GSY1135/GSY1135_Chr10.fasta ./Gen007.sort.bam > file.mpileup   

8.	Make depthchart files using pileup2dc
------
	cd src/pileup2dc

	gcc -o pileup2dc main.c

	mv pileup2dc ../../bin/

	PATH=$PATH:/home/fzhang/Research/rvd2/bin

	pileup2dc

9.	Estimate models for case + control
-----
10.	Test for variants between case + control -> vcf files
-----

11. http://www.biostars.org/p/19446/ 
----
	awk 'NR%2==1 { print $0 "/1" } ; NR%2==0 { print substr($0,0,length($0)/2) }' ../BARCODE/Gen007_unmatched.fastq > gen007_test1.fastq

	awk 'NR%2==1 { print $0 "/2" } ; NR%2==0 { print substr($0,length($0)/2) }' ../BARCODE/Gen007_unmatched.fastq > gen007_test2.fastq



------------------
###Notes1: http://biobits.org/samtools_primer.html
-------------
###Notes2: control + A + D
------------
		  screen -S/ -ls/ -r
		  ls -lh 



Prepared by: _________ _Fan Zhang______ Date: ____________12/20/2013____________


Witnessed by: ________________________