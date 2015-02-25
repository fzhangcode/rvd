#!/bin/sh
FASTAFILE = /research/pjflaherty/flahertylab/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
BAMDIR = /research/pjflaherty/flahertylab/freeze/HCC1187/
REGION = chr7:154735400-154794682
PILEUP2DC = ../../bin/pileup2dc
RVD2 = ../../src/rvd27/rvd27.py

echo 'pileup'
samtools mpileup -C 50 -d 100000  -r $(REGION) -f $(FASTAFILE) $(BAMDIR)/HCC1187BL_S1.bam > HCC1187BL_S1.pileup
samtools mpileup -C 50 -d 100000  -r $(REGION) -f $(FASTAFILE) $(BAMDIR)/HCC1187CL_S1.bam > HCC1187C_S1.pileup

echo 'pileup2dc'
$(PILEUP2DC) HCC1187BL_S1.pileup > HCC1187BL_S1.dc
$(PILEUP2DC) HCC1187C_S1.pileup > HCC1187C_S1.dc

echo 'Running RVD2 in the HCC1187 subject PAXIP1 gene dataset.'
echo 'Please provide the directory to the control and case hdf5 files.'

echo 'Generating control and case hdf5 file from depth chart files:'
python $(RVD2) gibbs -o HCC1187BL_S1 HCC1187BL_S1.dc -p 20 -s 19860522
python $(RVD2) gibbs -o HCC1187C_S1 HCC1187C_S1.dc -p 20 -s 19860522

echo 'Performing hypothesis test (germline test and somatic test)'
python $(RVD2) germline_test HCC1187BL_S1.hdf5 
python $(RVD2) somatic_test HCC1187BL_S1.hdf5 HCC1187C_S1.hdf5
echo 'Done.'
