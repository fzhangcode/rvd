#!/bin/sh

WORK_DIR=$(dirname $0)
STRELKA_DIR=../../bin/strelka_workflow
DATA_DIR=../2013-08-06_Downsample_Read_Depth/bam/10
FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

# Copy configuration ini file to a local copy
echo Copy configuration ini file to local
cp $STRELKA_DIR/etc/strelka_config_eland_default.ini config.ini

OUTPUT_DIR=./myAnalysis

samtools index $DATA_DIR/20100916_c3_p1.07_CGT.sorted.bam
samtools index $DATA_DIR/20100916_c1_p1.02_ACT.sorted.bam

# Configure
$STRELKA_DIR/bin/configureStrelkaWorkflow.pl \
	--normal=$DATA_DIR/20100916_c1_p1.02_ACT.sorted.bam \
	--tumor=$DATA_DIR/20100916_c3_p1.07_CGT.sorted.bam \
	--ref=$FASTAFILE \
	--config=config.ini \
	--output-dir=./myAnalysis
