#!/bin/sh
BAMDIR=../../data/Synthetic_BAM_files
BAMFILE=$BAMDIR/*.bam

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
VARSCAN2=../../bin/VarScan.v2.3.4.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)
MAXDEPTH=100000

for j in ${J[@]:0:4}
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

	VCFDIR=vcf/$DFRAC
	mkdir -p $VCFDIR
	PILEUPDIR=pileup/$DFRAC
	mkdir -p $PILEUPDIR

	SORT=($SORTBAM0_1 $SORTBAM0_3 $SORTBAM1_0 $SORTBAM10_0 $SORTBAM100_0)
	VCF=(vcf0_1 vcf0_3 vcf1_0 vcf10_0 vcf100_0)
	PILEUP=(control_case0_1 control_case0_3 control_case1_0 control_case10_0 control_case100_0)

	N=(0 1 2 3 4)
	for i in ${N[@]:0:5}
	do	
		let h=6\*i
		echo -------------------------------------------------------
		echo $DFRAC $PILEUP[$i]
		SORTBAM=($SORTCONTROL ${SORT[@]:$h:6})
		if [ -f $PILEUPDIR/${PILEUP[$i]}.pileup ]
			then
				echo File  $PILEUPDIR/${PILEUP[$i]}.pileup exists already 
			else
				samtools mpileup -C 50 -d 1000000 -f $FASTAFILE  ${SORTBAM[@]:0:12} > $PILEUPDIR/${PILEUP[$i]}.pileup
          	 fi
		
		echo -------------------------------------------------------
		echo SNP calling using VarScan2_Somatic
		if [ -f $VCFDIR/${VCF[$i]}.snp ]
			then
				echo File  $VCFDIR/${VCF[$i]}.snp exists already
			else
				java -jar $VARSCAN2 somatic $PILEUPDIR/${PILEUP[$i]}.pileup $VCFDIR/${VCF[$i]} --mpileup 1
					#--normal-purity 1.00
					#--tumor-purity 1.00
		fi 
	done
done
