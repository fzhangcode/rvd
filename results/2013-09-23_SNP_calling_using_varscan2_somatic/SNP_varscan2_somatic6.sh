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

	DOWNDIR=../2013-08-06_Downsample_Read_Depth/bam/$DFRAC
	
	VCFDIR=vcf6/$DFRAC
	mkdir -p $VCFDIR
	PILEUPDIR=pileup6/$DFRAC
	mkdir -p $PILEUPDIR


	controlbam=$DOWNDIR/20100916_c1_p1.02_ACT.sorted.bam
	case100_0bam=$DOWNDIR/20100916_c3_p1.07_CGT.sorted.bam
	echo -------------------------------------------------------
	echo pileup 
	if [ -f $PILEUPDIR/control_case100_0.pileup ]
		then
			echo File $PILEUPDIR/control_case100_0.pileup exists already
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $controlbam $case100_0bam > $PILEUPDIR/control_case100_0.pileup 
	fi


	echo -------------------------------------------------------
	echo call VarScan somatic		
	java -jar $VARSCAN2 somatic $PILEUPDIR/control_case100_0.pileup $VCFDIR/control_case100_0 \
		--mpileup 1 \
		--min-var-freq 0.00001 \
		--normal-purity 1.00 \
		--tumor-purity 1.00

done

