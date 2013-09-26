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
	
	SORTCONTROL=$(ls $DOWNDIR/20100916_c1_p?.02*.sorted.bam $DOWNDIR/20100916_c1_p?.04*.sorted.bam $DOWNDIR/20100916_c1_p?.05*.sorted.bam)
	SORTBAM0_1=$(ls $DOWNDIR/20100916_c1_p?.07*.sorted.bam  $DOWNDIR/20100916_c1_p?.12*.sorted.bam $DOWNDIR/20100916_c1_p?.14*.sorted.bam)
	SORTBAM0_3=$(ls $DOWNDIR/20100916_c2_p?.02*.sorted.bam  $DOWNDIR/20100916_c2_p?.04*.sorted.bam $DOWNDIR/20100916_c2_p?.05*.sorted.bam)
	SORTBAM1_0=$(ls $DOWNDIR/20100916_c2_p?.07*.sorted.bam  $DOWNDIR/20100916_c2_p?.12*.sorted.bam $DOWNDIR/20100916_c2_p?.14*.sorted.bam)
	SORTBAM10_0=$(ls $DOWNDIR/20100916_c3_p?.02*.sorted.bam  $DOWNDIR/20100916_c3_p?.04*.sorted.bam $DOWNDIR/20100916_c3_p?.05*.sorted.bam)
	SORTBAM100_0=$(ls $DOWNDIR/20100916_c3_p?.07*.sorted.bam  $DOWNDIR/20100916_c3_p?.12*.sorted.bam $DOWNDIR/20100916_c3_p?.14*sorted.bam)

	VCFDIR=vcf1/$DFRAC
	mkdir -p $VCFDIR
	PILEUPDIR=pileup1/$DFRAC
	mkdir -p $PILEUPDIR
		
	echo -------------------------------------------------------
	echo Pileup 6 control bam files
	if [ -f $PILEUPDIR/control.pileup ]
		then
			echo $PILEUPDIR/control.pileup exists
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $SORTCONTROL > $PILEUPDIR/control.pileup 
	fi

	echo -------------------------------------------------------
	echo Pileup 6 case100_0 bam files

	if [ -f $PILEUPDIR/case100_0.pileup ]
		then
			echo $PILEUPDIR/case100_0.pileup exists
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $SORTBAM100_0> $PILEUPDIR/case100_0.pileup 
	fi

	
	echo -------------------------------------------------------
	echo call VarScan somatic	
	java -jar $VARSCAN2 somatic  $PILEUPDIR/control.pileup $PILEUPDIR/case100_0.pileup $VCFDIR/control_case100_0 \
		--min-var-freq 0.00001 \
		--normal-purity 1.00 \
		--tumor-purity 1.00

done

