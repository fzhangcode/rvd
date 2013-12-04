#!/bin/sh

SCRIPT_DIR=$(dirname $0)

STRELKA_DIR=../../bin/strelka_workflow

FASTAFILE=../../data/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

DEFAULTCONFIG=$STRELKA_DIR/etc/strelka_config_eland_default.ini
BAMCONTROL=../../data/HCC1187/HCC1187BL_S1.bam
BAMTUMOR=../../data/HCC1187/HCC1187C_S1.bam

REGION=chr7:154738059-154782774

SUBBAMCONTROL=subHCC1187BL_S1.bam
SUBBAMTUMOR=subHCC1187C_S1.bam

# Use samtools to subtract HCC1187 PAXIP1 gene
samtools view -h -b $BAMCONTROL $REGION > $SUBBAMCONTROL
samtools index $SUBBAMCONTROL
samtools view -h -b $BAMTUMOR $REGION > $SUBBAMTUMOR
samtools index $SUBBAMTUMOR
echo done

echo ==========
WORKDIR=work
mkdir -p $WORKDIR
echo -------------------------------------------------------
echo Copy configuration ini file to work directory
cp $DEFAULTCONFIG $WORKDIR/config.ini
echo =======
echo SNP calling using strelka
		

echo ====
echo ../../data/HCC1187/*.bam

# Configure
$STRELKA_DIR/bin/configureStrelkaWorkflow.pl \
	--normal=$SUBBAMCONTROL \
	--tumor=$SUBBAMTUMOR \
	--ref=$FASTAFILE \
	--config=$WORKDIR/config.ini \
	--output-dir=$WORKDIR/myAnalysis
echo done
make -C $WORKDIR/myAnalysis

