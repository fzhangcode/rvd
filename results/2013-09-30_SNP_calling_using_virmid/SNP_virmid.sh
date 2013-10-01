#!/bin/sh
BAMDIR=../../data/Synthetic_BAM_files
BAMFILE=$BAMDIR/*.bam

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
VIRMID=../2013-09-30_SNP_calling_using_virmid/virmid-1.0.2/Virmid.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)

for j in ${J[@]:0:1}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	echo -------------------------------------------------------

	DOWNDIR=../2013-08-06_Downsample_Read_Depth/bam/$DFRAC
	
	VCFDIR=vcf/$DFRAC
	mkdir -p $VCFDIR


	controlbam=$DOWNDIR/20100916_c1_p1.04*.sorted.bam
	case100_0bam=$DOWNDIR/20100916_c3_p1.12*.sorted.bam

	echo -------------------------------------------------------
	echo Index the sorted bam files
	samtools index $controlbam
	samtools index $case100_0bam
	samtools faidx plasmid.fa
	echo -------------------------------------------------------
	echo call Virmid	
	java -jar $VIRMID -R plasmid.fa -D $case100_0bam -N $controlbam -w $VCFDIR/ -f

done

