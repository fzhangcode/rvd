#!/bin/sh
FASTAFILE=../../../../../freeze/HCC1187/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
BAMCONTROL=../../../../../freeze/HCC1187/HCC1187BL_S1.bam
BAMCASE=../../../../../freeze/HCC1187/HCC1187C_S1.bam
VARSCAN2=VarScan.v2.3.4.jar

echo -------------------------
echo extract PAXIP1 section of bam files
samtools view -u $BAMCONTROL chr7:154735400-154794682 > HCC1187BL_S1.part.bam
samtools view -u $BAMCASE chr7:154735400-154794682 > HCC1187C_S1.part.bam

echo -------------------------
echo index PAXIP1 section of bam files
samtools index HCC1187BL_S1.part.bam
samtools index HCC1187C_S1.part.bam

echo -------------------------
echo sort bam files
samtools sort HCC1187BL_S1.part.bam HCC1187BL_S1.part.sorted
samtools sort HCC1187C_S1.part.bam HCC1187C_S1.part.sorted

echo -------------------------
echo index sorted PAXIP1 section of bam files
samtools index HCC1187BL_S1.part.sorted.bam
samtools index HCC1187C_S1.part.sorted.bam


echo -------------------------
echo mpileup
samtools mpileup -C50 -f $FASTAFILE -r chr7:154735400-154794682 HCC1187BL_S1.part.sorted.bam > HCC1187BL_S1.part.pileup 
samtools mpileup -C50 -f $FASTAFILE -r chr7:154735400-154794682 HCC1187C_S1.part.sorted.bam > HCC1187C_S1.part.pileup

echo -------------------------
echo SNP calling using VARSCAN2 somatic
java -jar $VARSCAN2 somatic HCC1187BL_S1.part.pileup HCC1187C_S1.part.pileup ./HCC1187 \
					--min-var-freq 0.00001\
					--normal-purity 1.00 \
					--tumor-purity 1.00
