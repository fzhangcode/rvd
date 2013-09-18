#!/bin/sh
BAMDIR=../../data/Synthetic_BAM_files
BAMFILE=$BAMDIR/*.bam

PICARD=/usr/local/pjf/picard-tools-1.96
FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

echo $BAM0_1
echo $BAM0_3
echo $BAM1_0
echo $BAM10_0
echo $BAM100_0

DRATE=1
DFRAC=100000

echo -------------------------------------------------------------------
echo Start Downsampling and sorting


DOWNDIR=downsample/$DFRAC
mkdir -p $DOWNDIR
for f in $BAMFILE
do
	filename=${f##*/}
	Input=$f
	DownOutput=${DOWNDIR%%/}/$filename
	SortOutput=${DownOutput//bam/sorted}
	echo $Output

	if [ -f $DownOutput ]
		then 
			echo File $filename exists already
		else
			echo Downsampling $filename
			java -Xmx2g -jar $PICARD/DownsampleSam.jar INPUT=$Input OUTPUT=$DownOutput RANDOM_SEED=null PROBABILITY=$DRATE VALIDATION_STRINGENCY=SILENT
			samtools sort $Output
	fi
	
	if [ -f $SortOutput.bam ]
		then
			echo File $SortOutput.bam exists already
		else
			samtools sort $DownOutput $SortOutput
	fi
			
done

SORTBAM0_1=$(ls $DOWNDIR/20100916_c1_p?.07*.sorted.bam  $DOWNDIR/20100916_c1_p?.12*.sorted.bam $DOWNDIR/20100916_c1_p?.14*.sorted.bam)
SORTBAM0_3=$(ls $DOWNDIR/20100916_c2_p?.02*.sorted.bam  $DOWNDIR/20100916_c2_p?.04*.sorted.bam $DOWNDIR/20100916_c2_p?.05*.sorted.bam)
SORTBAM1_0=$(ls $DOWNDIR/20100916_c2_p?.07*.sorted.bam  $DOWNDIR/20100916_c2_p?.12*.sorted.bam $DOWNDIR/20100916_c2_p?.14*.sorted.bam)
SORTBAM10_0=$(ls $DOWNDIR/20100916_c3_p?.02*.sorted.bam  $DOWNDIR/20100916_c3_p?.04*.sorted.bam $DOWNDIR/20100916_c3_p?.05*.sorted.bam)
SORTBAM100_0=$(ls $DOWNDIR/20100916_c3_p?.07*.sorted.bam  $DOWNDIR/20100916_c3_p?.12*.sorted.bam $DOWNDIR/20100916_c3_p?.14*sorted.bam)

echo $SORTBAM0_1

VCFDIR=vcf/$DFRAC
mkdir -p $VCFDIR
SORT=($SORTBAM0_1 $SORTBAM0_3 $SORTBAM1_0 $SORTBAM10_0 $SORTBAM100_0)
VCF=(vcf0_1 vcf0_3 vcf1_0 vcf10_0 vcf100_0)
echo ---------------------------------------
N=(0 1 2 3 4)
echo ${N[@]:0:5}
for i in ${N[@]:0:5}
do
	let h=6\*i
	echo ${SORT[@]:$h:6}
	echo ---------------------------------------
	echo SNP Calling
	samtools mpileup -uf $FASTAFILE ${SORT[@]:$h:6} | bcftools view -bvcg - > $VCFDIR/${VCF[$i]}.bcf
	bcftools view VCFDIR/${VCF[$i]}.bcf > $VCFDIR/${VCF[$i]}.vcf 
done
