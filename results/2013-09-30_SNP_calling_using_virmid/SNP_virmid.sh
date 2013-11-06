#!/bin/sh

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
VIRMID=../2013-09-30_SNP_calling_using_virmid/virmid-1.0.2/Virmid.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)

for j in ${J[@]:0:1}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	DOWNDIR=../2013-10-02_problematic_header_removal/bam/$DFRAC

	VCFDIR=vcf/$DFRAC
	mkdir -p $VCFDIR

	controlbam=$DOWNDIR/20100916_c1_p1.02_ACT.reheadered.sorted.bam
	case100_0bam=$DOWNDIR/20100916_c3_p1.07_CGT.reheadered.sorted.bam

	echo -------------------------------------------------------
	echo call Virmid	
	java -jar $VIRMID \
		-R $FASTAFILE \
		-D $case100_0bam \
		-N $controlbam \
		-r 400 \
		-q 0 \
		-w $VCFDIR/ \
		-v 3 \
		-af

done

