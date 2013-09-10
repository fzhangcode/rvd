#!/bin/sh

samtools mpileup -ug -d 100000 -f ../../data/Synthetic_BAM_files/plasmid.fa bam/20100916_c3_p1.05_CAT.sorted.bam bam/20100916_c3_p2.05_CAT.sorted.bam | bcftools view -bvcg - > var.raw.bcf

bcftools view var.raw.bcf | /home/pjflaherty/flahertylab/apps/samtools/bcftools/vcfutils.pl varFilter -D100000 > var.flt.vcf
