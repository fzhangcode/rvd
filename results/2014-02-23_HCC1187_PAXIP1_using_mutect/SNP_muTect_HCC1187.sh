#!/bin/sh

FASTAFILE=../../data/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

MUTECT=muTect.jar

BAMCONTROL=../../data/HCC1187/HCC1187BL_S1.bam
BAMCASE=../../data/HCC1187/HCC1187C_S1.bam


WORKDIR=work/
mkdir -p $WORKDIR

java -Xmx2g -jar $MUTECT \
	--analysis_type MuTect \
	--reference_sequence $FASTAFILE \
	--intervals chr7:154738059-154782774 \
	--input_file:normal $BAMCONTROL \
	--input_file:tumor $BAMCASE \
	--out $WORKDIR/HCC.call_stats.txt \
	--coverage_file $WORKDIR/HCC.coverage.wig.txt
