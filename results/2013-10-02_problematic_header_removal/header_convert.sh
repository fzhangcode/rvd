#!/bin/sh

FASTAFILE=../../data/Synthetic_BAM_files/plasmid.fa

DRATELIST=(0.1 0.01 0.001 0.0001)
DILUTIONLIST=(0.1 0.3 1.0 10.0 100.0)
J=(0 1 2 3)

for j in ${J[@]:0:4}
do	
	DRATE=${DRATELIST[$j]}
	DFRAC=$(echo $DRATE| awk '{printf "%2.0f\n",1/$1}')

	DOWNDIR=../2013-08-06_Downsample_Read_Depth/bam/$DFRAC
	
	OUTPUTDIR=bam/$DFRAC
	mkdir -p $OUTPUTDIR
	
	BAMFILE=$DOWNDIR/*.sorted.bam	
	for f in $BAMFILE
		do
			filename=${f##*/}
			basename=${filename%.sorted.bam}
			output=$OUTPUTDIR/$basename.reheadered.sorted.bam
			echo ----------------------------------
			echo Rehearder $f
			if [ -f $output ]
			then
				echo $output exists already
			else
				samtools view -H $f| \
					sed -e '3d'| \
					samtools reheader - $f > \
					$output
			fi
			
			echo ----------------------------------
			echo Index $output

			if [ -f $output.bai ]
                        then    
                                echo $output.bai exists already
                        else
				samtools index $output
                        fi
		done
done
