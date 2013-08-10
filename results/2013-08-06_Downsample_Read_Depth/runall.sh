#!/bin/sh
echo "Run program to downsample synthetic bam files using Picard"
echo "Please input desired downsampling rate"
read dRate
java -Xmx2g -jar \picard-tools-1.96\DownsampleSam.jar I=../../data/Synthetic_BAM_files/20100916_c1_p1.02_ACT.bam O=../../data/Synthetic/20100916_c1_p1.02_ACT_ds.bam R=100 P=dRate

