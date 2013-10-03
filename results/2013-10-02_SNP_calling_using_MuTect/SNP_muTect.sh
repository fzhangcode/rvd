#!/bin/sh

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

DOWNDIR=../2013-09-13_SNP_calling_using_GATK/format/10
controlbam=$DOWNDIR/20100916_c1_p1.02_ACT.sorted.fixed.bam
case100_0bam=$DOWNDIR/20100916_c3_p1.07_CGT.sorted.fixed.bam

java -Xmx2g -jar muTect-1.1.4.jar \
	--analysis_type MuTect \
	--reference_sequence $FASTAFILE \
	--input_file:normal $controlbam \
	--input_file:tumor $case100_0bam \
	--out example.call_stats.txt \
	--coverage_file example.coverage.wig.txt \
