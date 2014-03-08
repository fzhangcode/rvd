2014-02-24 SNP calling using RVD in synthetic dataset
==============================

Purpose
------------
Compare variant calling performance of RVD and RVD2 in synthetic dataset.

Conclusions
-----------------
RVD is able to call variants with synthetic data at 10X downsampling rate. The program failed when the downsampling rate is 100X or higher. 

Background
----------------

Materials and Equipment
------------------------------
The RVD program was downloaded from [RVD program files](http://hamachi.stanford.edu/publication-material/rvd/matlab/Program_Files.zip "RVD").

Clinical bam files are available in 
`/flahertylab/freeze/HCC1187`; reference fasta file is in `/flahertylab/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa`

More information regarding RVD can be found in [Algorithm User Guide for RVD](http://www.webcitation.org/query.php?url=http://dna-discovery.stanford.edu/software/rvd/&refdoi=10.1186/1756-0500-6-206)

Experimental Protocol
---------------------------
Please run matlab script RVD.m for the overall process. 


Results
-----------
Variants were called and stored in `output` folder for the 10x downsampling. The results match closely with those reported in the NAR paper.

The RVD command line program errored out for downsampling rates at 100x or higher due to insufficient depth in the error bases.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
