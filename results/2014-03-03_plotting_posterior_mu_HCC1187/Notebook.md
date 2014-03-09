2014-02-23 Histogram of empirical mu of positions called by MuTect in HCC1187 data
==============================

Purpose
------------
Plot the histogram of of empirical mu of positions called by MuTect in HCC1187 data to find out how mucase and mucontrol distributed in positions RVD2 missed but MuTect called.

Conclusions
-----------------
 

Background
----------------
From experiment `./2014-02-23_HCC1187_PAXIP1_using_mutect`, we noticed that RVD2 missed some positions Mutect called in clinical sample HCC1187 PAXIP1 gene. The depth matrix for these positions are shown in below. 

![]('HCC_call_stats.png')

Fig.1  Depth matrix for positions called by Mutect in HCC1187.

In order to figure out why these positions are missed, the empirical posterior distributions of mu _j (histogram) are plotted for visulization.


Materials and Equipment
------------------------------
The positions in Figure 1. can be generated in experiment `2014-02-23_HCC1187_PAXIP1_using_mutect`.

Please run the function `muplot.py` to generate the histograms. The control and case files are generated in folder `./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome`.


[How to run MuTect](http://www.broadinstitute.org/cancer/cga/mutect_run) instruction.

Clinical bam files are available in 
`/flahertylab/freeze/HCC1187`; reference fasta file is in `/flahertylab/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa`


Experimental Protocol
---------------------------
Please run the function `muplot.py` to generate the histograms.

Positions that are not called by RVD2 but called by Mutect:
![]('position154749704.png')
![]('position154766700.png')
![]('position154766732.png')
![]('position154777014.png')
![]('position154777118.png')

Positions that are called by both:
![]('position154743899.png')
![]('position154753635.png')
![]('position154754371.png')
![]('position154758813.png')
![]('position154760439.png')
![]('position154780960.png')

Results
-----------

The histogram of mu for the called positions are 



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
