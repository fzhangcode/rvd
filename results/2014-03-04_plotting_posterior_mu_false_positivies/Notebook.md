2014-03-04 Plot the posterior distribution for some typical false positive positions
==============================

Purpose
------------
Plot the histogram of of empirical mu of positions called by falsely called by RVD2 in HCC1187 dataset to actually see how the posterior looks like and why the false calls happen.

Conclusions
-----------------
It can be seen that there is problem in the sampling procedure which leads to the false calls. There is still no efficient solution to this problem.
 

Background
----------------
The algorithm has a limitation in sampling mu when the non-reference r=0 or r=n, which also means when mu is nearly 0 or 1. Actually, without chi2 test, there are positions called with depth chart shown in below:


![]('depthchart_for_false_calls.png')

Fig.1  Depth matrix for positions falsely called by RVD2 in HCC1187.

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
The positions in Figure 1. can be generated in experiment `2014-02-23_HCC1187_PAXIP1_using_mutect`.

Please run the function `muplot.py` to generate the histograms. The control and case files are generated in folder `./../2013-12-20_experiment_set_gibbs_Qsd_mu_1_mu_over_10_minus_mu0/HCC1187_PAXIP1_genome`




Results
-----------
Please see the histogram for this positions. 
![]('position154782559.png')
![]('position154782666.png')
![]('position154782703.png')
![]('position154782715.png')
![]('position154782742.png')

Figure 2. Histogram for the posterior mu for representative falsely called positions. The left panel shows histograms for mucase and mucontrol independently. The right panel shows histogram for mucase-mucontrol.

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
