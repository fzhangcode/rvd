#!/bin/sh
FASTAFILE=/research/pjflaherty/flahertylab/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
MUTECT=/research/pjflaherty/flahertylab/fzhang/Research/rvd/bin/muTect-1.1.4.jar
BAMCONTROL=/research/pjflaherty/flahertylab/freeze/HCC1187/HCC1187BL_S1.bam
BAMCASE=/research/pjflaherty/flahertylab/freeze/HCC1187/HCC1187C_S1.bam
WORKDIR=work_muTect/
#mkdir -p $WORKDIR
java -Xmx2g -jar ${MUTECT} \
	--analysis_type MuTect \
	--reference_sequence $FASTAFILE \
	--intervals chr7:154735400-154794682 \
	--input_file:normal $BAMCONTROL \
	--input_file:tumor $BAMCASE \
	--out $WORKDIR/HCC.call_stats.txt \
	--coverage_file $WORKDIR/HCC.coverage.wig.txt
