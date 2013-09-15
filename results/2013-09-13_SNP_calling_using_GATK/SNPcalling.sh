#!/bin/sh

GATK=/home/pjflaherty/flahertylab/apps/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar

PICARDFORMAT=/usr/local/pjf/picard-tools-1.96/AddOrReplaceReadGroups.jar

BAMDIR=../../data/Synthetic_BAM_files/20100916_c2_p1.07_CGT.bam

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

SORTBAMDIR=./sort

FILE=20100916_c2_p1.07_CGT.bam
SORTEDFILE=${FILE//bam/sorted}
echo $SORTEDFILE
# sort the BAM files
if [ -f 20100916_c2_p1.07_CGT.sorted.bam ]
then
	echo the file exists
else	
	samtools sort $BAMDIR $SORTEDFILE
fi

# format read group using picard
if [ -f 20100916_c2_p1.07_CGT.sorted.fixed.bam ]
then
	echo the file exists
else
	java -jar $PICARDFORMAT I=20100916_c2_p1.07_CGT.sorted.bam O=20100916_c2_p1.07_CGT.sorted.fixed.bam SORT_ORDER=coordinate RGID=H1N1 RGLB=H1N1 RGPL=illumina RGPU=H1N1 RGSM=H1N1 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT
fi

# find realignment regions
if [ -f 20100916_c2_p1.07_CGT.realign.intervals ]
then
	echo the file exists
else
	java -Xmx1g -jar $GATK -T RealignerTargetCreator -R $FASTAFILE -I 20100916_c2_p1.07_CGT.sorted.fixed.bam -o 20100916_c2_p1.07_CGT.realign.intervals
fi

# realignment
if [ -f 20100916_c2_p1.07_CGT.realigned.bam ]
then
	echo the file exists
else

	java -Xmx4g -jar $GATK \
		-I 20100916_c2_p1.07_CGT.sorted.fixed.bam \
		-R $FASTAFILE \
		-T IndelRealigner \
	-targetIntervals 20100916_c2_p1.07_CGT.realign.intervals \
	-o 20100916_c2_p1.07_CGT.realigned.bam
fi

#SNP calling
java -jar $GATK \
	-R $FASTAFILE \
	-T UnifiedGenotyper \
	-I 20100916_c2_p1.07_CGT.realigned.bam \
	-o var_GATK.vcf \
	--output_mode EMIT_VARIANTS_ONLY \
