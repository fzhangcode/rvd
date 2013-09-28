#!/bin/sh
BAMDIR=../../data/Synthetic_BAM_files
BAMFILE=$BAMDIR/*.bam

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa
VARSCAN2=../../bin/VarScan.v2.3.4.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
DILUTIONLIST=(0.1 0.3 1.0 10.0 100.0)
J=(0 1 2 3)

for j in ${J[@]:0:4}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	DOWNDIR=../2013-08-06_Downsample_Read_Depth/bam/$DFRAC
	
	SORTCONTROL=$DOWNDIR/20100916_c1_p1.02*.sorted.bam
	SORTBAM0_1=$DOWNDIR/20100916_c1_p1.07*.sorted.bam
	SORTBAM0_3=$DOWNDIR/20100916_c2_p1.02*.sorted.bam
	SORTBAM1_0=$DOWNDIR/20100916_c2_p1.07*.sorted.bam
	SORTBAM10_0=$DOWNDIR/20100916_c3_p1.02*.sorted.bam
	SORTBAM100_0=$DOWNDIR/20100916_c3_p1.07*.sorted.bam

	SNPDIR=snp/$DFRAC
	mkdir -p $SNPDIR
	PILEUPDIR=pileup/$DFRAC
	mkdir -p $PILEUPDIR

	SORT=($SORTBAM0_1 $SORTBAM0_3 $SORTBAM1_0 $SORTBAM10_0 $SORTBAM100_0)
	SNP=(snp0_1 snp0_3 snp1_0 snp10_0 snp100_0)
	PILEUP=(case0_1 case0_3 case1_0 case10_0 case100_0)

	echo -------------------------------------------------------
	echo pileup $SORTCONTROL 
	if [ -f $PILEUPDIR/control.pileup ]
		then
			echo File $PILEUPDIR/control.pileup exists already
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $SORTCONTROL > $PILEUPDIR/control.pileup 
	fi


	N=(0 1 2 3 4)
	for i in ${N[@]:0:5}
	do	
		let h=1\*i

		echo -------------------------------------------------------
		echo Pileup $PILEUPDIR/${PILEUP[$i]}
		if [ -f $PILEUPDIR/${PILEUP[$i]}.pileup ]
		then
			echo File $PILEUPDIR/${PILEUP[$i]}.pileup exists already
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE ${SORT[@]:$h:1} > $PILEUPDIR/${PILEUP[$i]}.pileup
		fi


		echo -------------------------------------------------------
		echo SNP calling using VarScan2_Somatic
		if [ -f $SNPDIR/${SNP[$i]}.snp ]
			then
				echo File  $SNPDIR/${SNP[$i]}.snp exists already
			else
				java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/${PILEUP[$i]}.pileup $SNPDIR/${SNP[$i]} \
					--min-var-freq 0.00001\
					--normal-purity 1.00 \
					--tumor-purity ${DILUTIONLIST[$i]} 
		fi 
	done
done
