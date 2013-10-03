#!/bin/sh

SCRIPT_DIR=$(dirname $0)

STRELKA_DIR=../../bin/strelka_workflow

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

DEFAULTCONFIG=$STRELKA_DIR/etc/strelka_config_eland_default.ini


DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)

for j in ${J[@]:0:4}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	DOWNDIR=../2013-10-02_problematic_header_removal/bam/$DFRAC
	
	BAMCONTROL=$DOWNDIR/20100916_c1_p1.02_ACT.reheadered.sorted.bam
	echo $BAMCONTROL
	BAM0_1=$DOWNDIR/20100916_c1_p1.07_CGT.reheadered.sorted.bam
	BAM0_3=$DOWNDIR/20100916_c2_p1.02_ACT.reheadered.sorted.bam
	BAM1_0=$DOWNDIR/20100916_c2_p1.07_CGT.reheadered.sorted.bam
	BAM10_0=$DOWNDIR/20100916_c3_p1.02_ACT.reheadered.sorted.bam
	BAM100_0=$DOWNDIR/20100916_c3_p1.07_CGT.reheadered.sorted.bam
	
	BAM=($BAM0_1 $BAM0_3 $BAM1_0 $BAM10_0 $BAM100_0)
	DILUTIONLIST=(0_1 0_3 1_0 10_0 100_0)

	echo -------------------------------------------------------

	N=(0 1 2 3 4)
	for i in ${N[@]:0:5}
	do	
		let h=1\*i
		WORKDIR=work/$DFRAC/${DILUTIONLIST[$i]}
		mkdir -p $WORKDIR
		echo -------------------------------------------------------
		echo Copy configuration ini file to work directory
		cp $DEFAULTCONFIG $WORKDIR/config.ini

		echo SNP calling using strelka
		echo $BAMCONTROL
		echo ${BAM[@]:$h:1}
		
		# Configure
		if [ -d $WORKDIR/myAnalysis ]
			then
				echo Analysis completed.
			else
				$STRELKA_DIR/bin/configureStrelkaWorkflow.pl \
					--normal=$BAMCONTROL \
					--tumor=${BAM[@]:$h:1} \
					--ref=$FASTAFILE \
					--config=$WORKDIR/config.ini \
					--output-dir=$WORKDIR/myAnalysis
				make -C $WORKDIR/myAnalysis
		fi
		echo $PWD

		cp $PWD/work/$DFRAC/${DILUTIONLIST[$i]}/myAnalysis/results/all.somatic.snvs.vcf $PWD/work/$DFRAC/vcf${DILUTIONLIST[$i]}.vcf
	done
done
