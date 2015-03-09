#!/bin/sh
BAMDIR=../2015-03-03_run_GATK_SS_FDR_downsample_depth_100/bam/
BAMFILE=$BAMDIR/*.bam

FASTAFILE=../2015-03-03_run_GATK_SS_FDR_downsample_depth_100/bam/plasmid.fa
VARSCAN2=../../bin/VarScan.v2.3.4.jar

#DRATELIST=(0.1 0.01 0.001 0.0001)
J=(0 1 2 3)

DRATELIST=(0.0001)
J=(0)

for j in ${J[@]:0:1}
do
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	echo -------------------------------------------------------

	DOWNDIR=../2015-03-03_run_GATK_SS_FDR_downsample_depth_100/bam/$DFRAC
	
	SORTCONTROL=$(ls $DOWNDIR/20100916_c1_p?.02*.sorted.bam $DOWNDIR/20100916_c1_p?.04*.sorted.bam $DOWNDIR/20100916_c1_p?.05*.sorted.bam)
	SORTBAM0_1=$(ls $DOWNDIR/20100916_c1_p?.07*.sorted.bam  $DOWNDIR/20100916_c1_p?.12*.sorted.bam $DOWNDIR/20100916_c1_p?.14*.sorted.bam)
	SORTBAM0_3=$(ls $DOWNDIR/20100916_c2_p?.02*.sorted.bam  $DOWNDIR/20100916_c2_p?.04*.sorted.bam $DOWNDIR/20100916_c2_p?.05*.sorted.bam)
	SORTBAM1_0=$(ls $DOWNDIR/20100916_c2_p?.07*.sorted.bam  $DOWNDIR/20100916_c2_p?.12*.sorted.bam $DOWNDIR/20100916_c2_p?.14*.sorted.bam)
	SORTBAM10_0=$(ls $DOWNDIR/20100916_c3_p?.02*.sorted.bam  $DOWNDIR/20100916_c3_p?.04*.sorted.bam $DOWNDIR/20100916_c3_p?.05*.sorted.bam)
	SORTBAM100_0=$(ls $DOWNDIR/20100916_c3_p?.07*.sorted.bam  $DOWNDIR/20100916_c3_p?.12*.sorted.bam $DOWNDIR/20100916_c3_p?.14*sorted.bam)

	VCFDIR=snp3/$DFRAC
	mkdir -p $VCFDIR
	PILEUPDIR=pileup3/$DFRAC
	mkdir -p $PILEUPDIR


        echo ------------100-------------------------------------------
        echo merge control bam files
        if [ -f $PILEUPDIR/control.sorted.bam ]
                then
                        echo File $PILEUPDIR/control.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/control.sorted.bam $SORTCONTROL
        fi
 
        echo -------------------------------------------------------
        echo merge case bam files
        if [ -f $PILEUPDIR/case100_0.sorted.bam ]
                then
                        echo File $PILEUPDIR/case100_0.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/case100_0.sorted.bam $SORTBAM100_0
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/control.sorted.bam
        if [ -f $PILEUPDIR/control.pileup ]
                then
                        echo File $PILEUPDIR/control.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam  > $PILEUPDIR/control.pileup
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/case100_0.bam 
        if [ -f $PILEUPDIR/case100_0.pileup ]
                then
                        echo File $PILEUPDIR/case100_0.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case100_0.sorted.bam  > $PILEUPDIR/case100_0.pileup  
        fi

        java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case100_0.pileup $VCFDIR/control_case100_0 \
                --min-var-freq 0.00001 \
                --normal-purity 1.00 \
                --tumor-purity 1.00



	echo ----------------10---------------------------------------
	echo merge control bam files
	if [ -f $PILEUPDIR/control.sorted.bam ]
		then
			echo File $PILEUPDIR/control.sorted.bam exists already
		else
			samtools merge -u $PILEUPDIR/control.sorted.bam $SORTCONTROL
	fi
 
	echo -------------------------------------------------------
	echo merge case bam files
	if [ -f $PILEUPDIR/case10_0.sorted.bam ]
		then
			echo File $PILEUPDIR/case10_0.sorted.bam exists already
		else
			samtools merge -u $PILEUPDIR/case10_0.sorted.bam $SORTBAM10_0
	fi

	echo -------------------------------------------------------
	echo pileup $PILEUPDIR/control.sorted.bam 
	if [ -f $PILEUPDIR/control.pileup ]
		then
			echo File $PILEUPDIR/control.pileup exists already
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam  > $PILEUPDIR/control.pileup 
	fi

	echo -------------------------------------------------------
	echo pileup $PILEUPDIR/case10_0.bam 
	if [ -f $PILEUPDIR/case10_0.pileup ]
		then
			echo File $PILEUPDIR/case10_0.pileup exists already
		else
			samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case10_0.sorted.bam  > $PILEUPDIR/case10_0.pileup 
	fi
		
	java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case10_0.pileup $VCFDIR/control_case10_0 \
		--min-var-freq 0.00001 \
		--normal-purity 1.00 \
		--tumor-purity 1.00

                
        echo ----------------1---------------------------------------
        echo merge control bam files
        if [ -f $PILEUPDIR/control.sorted.bam ]
                then
                        echo File $PILEUPDIR/control.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/control.sorted.bam $SORTCONTROL
        fi
 
        echo -------------------------------------------------------
        echo merge case bam files
        if [ -f $PILEUPDIR/case1_0.sorted.bam ]
                then
                        echo File $PILEUPDIR/case1_0.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/case1_0.sorted.bam $SORTBAM1_0
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/control.sorted.bam
        if [ -f $PILEUPDIR/control.pileup ]
                then
                        echo File $PILEUPDIR/control.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam  > $PILEUPDIR/control.pileup
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/case1_0.bam
        if [ -f $PILEUPDIR/case1_0.pileup ]
                then
                        echo File $PILEUPDIR/case1_0.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case1_0.sorted.bam  > $PILEUPDIR/case1_0.pileup
        fi

        java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case1_0.pileup $VCFDIR/control_case1_0 \
                --min-var-freq 0.00001 \
                --normal-purity 1.00 \
                --tumor-purity 1.00

        echo ---------------0_3---------------------------------------
        echo merge control bam files
        if [ -f $PILEUPDIR/control.sorted.bam ]
                then
                        echo File $PILEUPDIR/control.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/control.sorted.bam $SORTCONTROL
        fi
 
        echo -------------------------------------------------------
        echo merge case bam files
        if [ -f $PILEUPDIR/case0_3.sorted.bam ]
                then
                        echo File $PILEUPDIR/case0_3.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/case0_3.sorted.bam $SORTBAM0_3
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/control.sorted.bam
        if [ -f $PILEUPDIR/control.pileup ]
                then
                        echo File $PILEUPDIR/control.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam  > $PILEUPDIR/control.pileup
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/case0_3.bam
        if [ -f $PILEUPDIR/case0_3.pileup ]
                then
                        echo File $PILEUPDIR/case0_3.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case0_3.sorted.bam  > $PILEUPDIR/case0_3.pileup
        fi

        java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case0_3.pileup $VCFDIR/control_case0_3 \
                --min-var-freq 0.00001 \
                --normal-purity 1.00 \
                --tumor-purity 1.00

        echo ---------------0_1---------------------------------------
        echo merge control bam files
        if [ -f $PILEUPDIR/control.sorted.bam ]
                then
                        echo File $PILEUPDIR/control.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/control.sorted.bam $SORTCONTROL
        fi
 
        echo -------------------------------------------------------
        echo merge case bam files
        if [ -f $PILEUPDIR/case0_1.sorted.bam ]
                then
                        echo File $PILEUPDIR/case0_1.sorted.bam exists already
                else
                        samtools merge -u $PILEUPDIR/case0_1.sorted.bam $SORTBAM0_1
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/control.sorted.bam
        if [ -f $PILEUPDIR/control.pileup ]
                then
                        echo File $PILEUPDIR/control.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/control.sorted.bam  > $PILEUPDIR/control.pileup
        fi

        echo -------------------------------------------------------
        echo pileup $PILEUPDIR/case0_1.bam
        if [ -f $PILEUPDIR/case0_1.pileup ]
                then
                        echo File $PILEUPDIR/case0_1.pileup exists already
                else
                        samtools mpileup -C 50 -d 100000 -f $FASTAFILE $PILEUPDIR/case0_1.sorted.bam  > $PILEUPDIR/case0_1.pileup
        fi

        java -jar $VARSCAN2 somatic $PILEUPDIR/control.pileup $PILEUPDIR/case0_1.pileup $VCFDIR/control_case0_1 \
                --min-var-freq 0.00001 \
                --normal-purity 1.00 \
                --tumor-purity 1.00

done

