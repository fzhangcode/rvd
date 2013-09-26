#!/bin/sh
BAMDIR=../../data/Synthetic_BAM_files
BAMFILE=$BAMDIR/*.bam

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
VARSCAN2=../../bin/VarScan.v2.3.4.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)
MAXDEPTH=100000

for j in ${J[@]:0:1}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE $MAXDEPTH | awk '{printf "%2.0f\n",$1*$2}')

	echo -------------------------------------------------------

	
	
	DOWNDIR=../../data/Synthetic_BAM_files
	VCFDIR=vcf7/$DFRAC
	mkdir -p $VCFDIR
	PILEUPDIR=pileup7/$DFRAC
	mkdir -p $PILEUPDIR
	
	controlbam=$DOWNDIR/20100916_c1_p1.02_ACT.bam
	case100_0bam=$DOWNDIR/20100916_c3_p1.07_CGT.bam

	echo -------------------------------------------------------
	echo sort bam files
	samtools sort $controlbam $PILEUPDIR/control.sorted
	samtools sort $case100_0bam $PILEUPDIR/case100_0.sorted

	echo -------------------------------------------------------
	echo pileup 
	samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam > $PILEUPDIR/control.pileup 
	samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case100_0.sorted.bam > $PILEUPDIR/case100_0.pileup 

	echo -------------------------------------------------------
	echo call VarScan somatic		
	java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case100_0.pileup $VCFDIR/control_case100_0 \
		--min-var-freq 0.00001 \
		--normal-purity 1.00 \
		--tumor-purity 1.00

done

