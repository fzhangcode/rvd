#!/bin/sh

GATK=/home/pjflaherty/flahertylab/apps/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar
PICARD=/usr/local/pjf/picard-tools-1.96

BAMDIR=../../data/Synthetic_BAM_files
FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
BAMFILE=$BAMDIR/*.bam

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)
MAXDEPTH=100000

for j in ${J[@]:0:4}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE $MAXDEPTH | awk '{printf "%2.0f\n",$1*$2}')

	# Preprocessing
	echo Preprocessing----------------------

	BAMDIR=../2013-08-06_Downsample_Read_Depth/bam/$DFRAC
	BAMFILE=$BAMDIR/*.sorted.bam


	FORMATDIR=format/$DFRAC
	mkdir -p $FORMATDIR

	for f in $BAMFILE
	do
		filename=${f##*/}
		echo format read group using picard-------
		FormatOutput=$FORMATDIR/$filename
		FormatOutput=${FormatOutput//sorted.bam/sorted.fixed.bam}
		if [ -f $FormatOutput ]
			then
				echo File $FormatOutput exists already
			else
				java -jar $PICARD/AddOrReplaceReadGroups.jar \
					I=$f \
					O=$FormatOutput \
					SORT_ORDER=coordinate \
					RGID=H1N1 \
					RGLB=H1N1 \
					RGPL=illumina \
					RGPU=H1N1 \
					RGSM=H1N1 \
					CREATE_INDEX=True \
					VALIDATION_STRINGENCY=LENIENT
		fi
	
		echo find realignment regions----------------
		RealignInterval=${FormatOutput//sorted.fixed.bam/realign.intervals}
		if [ -f $RealignInterval ]
			then
				echo File $RealignInterval exists already
			else
				java -Xmx1g -jar $GATK \
					-T RealignerTargetCreator \
					-R $FASTAFILE \
					-I $FormatOutput \
					-o $RealignInterval
		fi

		echo realignment------------------------------
		RealignOutput=${FormatOutput//sorted.fixed.bam/realign.bam}
		if [ -f $RealignOutput ]
			then
				echo File $RealignOutput exists already
			else
				java -Xmx4g -jar $GATK \
					-I $FormatOutput \
					-R $FASTAFILE \
					-T IndelRealigner \
					-targetIntervals $RealignInterval \
					-o $RealignOutput
		fi
	done

	Realign0_1=$(ls $FORMATDIR/20100916_c1_p?.07*.realign.bam  $FORMATDIR/20100916_c1_p?.12*.realign.bam $FORMATDIR/20100916_c1_p?.14*.realign.bam)
	Realign0_3=$(ls $FORMATDIR/20100916_c2_p?.02*.realign.bam  $FORMATDIR/20100916_c2_p?.04*.realign.bam $FORMATDIR/20100916_c2_p?.05*.realign.bam)
	Realign1_0=$(ls $FORMATDIR/20100916_c2_p?.07*.realign.bam  $FORMATDIR/20100916_c2_p?.12*.realign.bam $FORMATDIR/20100916_c2_p?.14*.realign.bam)
	Realign10_0=$(ls $FORMATDIR/20100916_c3_p?.02*.realign.bam  $FORMATDIR/20100916_c3_p?.04*.realign.bam $FORMATDIR/20100916_c3_p?.05*.realign.bam)
	Realign100_0=$(ls $FORMATDIR/20100916_c3_p?.07*.realign.bam  $FORMATDIR/20100916_c3_p?.12*.realign.bam $FORMATDIR/20100916_c3_p?.14*.realign.bam)

	echo SNP Calling---------------------------
	VCFDIR=vcf/$DFRAC
	mkdir -p $VCFDIR
	REALIGN=($Realign0_1 $Realign0_3 $Realign1_0 $Realign10_0 $Realign100_0)
	VCF=(vcf0_1 vcf0_3 vcf1_0 vcf10_0 vcf100_0)
	echo ---------------------------------------
	N=(0 1 2 3 4)
	for i in ${N[@]:0:5}
	do
		if [ -f $VCFDIR/${VCF[$i]}.vcf ]
			then
				echo File  $VCFDIR/${VCF[$i]}.vcf exists already
			else
				let h=6\*i
				echo ---------------------------------------------
				java -jar $GATK \
				-R $FASTAFILE \
				-T UnifiedGenotyper \
				-I ${REALIGN[@]:$h:1} \
				-I ${REALIGN[@]:$h+1:1} \
				-I ${REALIGN[@]:$h+2:1} \
				-I ${REALIGN[@]:$h+3:1} \
				-I ${REALIGN[@]:$h+4:1} \
				-I ${REALIGN[@]:$h+5:1} \
				-o $VCFDIR/${VCF[$i]}.vcf \
				--output_mode EMIT_VARIANTS_ONLY
		fi
	done
done
