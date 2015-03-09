#!/bin/sh

FASTAFILE=../2015-03-03_run_GATK_SS_FDR_downsample_depth_100/bam/plasmid.fa

#DEFAULTCONFIG=$STRELKA_DIR/etc/strelka_config_eland_default.ini

MUTECT=../../bin/muTect-1.1.4.jar

DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)

for j in ${J[@]:0:4}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	#DOWNDIR=../2013-09-13_SNP_calling_using_GATK/format/$DFRAC
	DOWNDIR=../2015-03-03_run_GATK_SS_FDR_downsample_depth_100/format/$DFRAC
	
	BAMCONTROL=$DOWNDIR/20100916_c1_p2.04_ATT.sorted.fixed.bam
	BAM0_1=$DOWNDIR/20100916_c1_p2.07_CGT.sorted.fixed.bam
	BAM0_3=$DOWNDIR/20100916_c2_p2.02_ACT.sorted.fixed.bam
	BAM1_0=$DOWNDIR/20100916_c2_p1.12_GTT.sorted.fixed.bam
	BAM10_0=$DOWNDIR/20100916_c3_p1.04_ATT.sorted.fixed.bam
	BAM100_0=$DOWNDIR/20100916_c3_p1.12_GTT.sorted.fixed.bam
	
	BAM=($BAM0_1 $BAM0_3 $BAM1_0 $BAM10_0 $BAM100_0)
	DILUTIONLIST=(0_1 0_3 1_0 10_0 100_0)

	N=(0 1 2 3 4)
	for i in ${N[@]:0:5}
	do	
		let h=1\*i
		#WORKDIR=work/$DFRAC/${DILUTIONLIST[$i]}
		WORKDIR=work/$DFRAC/

		mkdir -p $WORKDIR
		echo -------------------------------------------------------
		
		# Call MuTect
		if [ -f $WORKDIR/${DILUTIONLIST[$i]}.call_stats.txt ]
			then
				echo Analysis completed.
			else
				java -Xmx2g -jar $MUTECT \
					--analysis_type MuTect \
					--reference_sequence $FASTAFILE \
					--input_file:normal $BAMCONTROL \
					--input_file:tumor ${BAM[@]:$h:1} \
					--out $WORKDIR/${DILUTIONLIST[$i]}.call_stats.txt \
					--coverage_file $WORKDIR/${DILUTIONLIST[$i]}.coverage.wig.txt
		fi
	done
done

echo -------------------------------------------------------
echo Convert snp files to vcf files...
python txt2vcf.py
