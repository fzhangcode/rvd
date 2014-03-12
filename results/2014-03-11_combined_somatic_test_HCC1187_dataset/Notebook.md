2014-03-11 combine somatic test in HCC1187 dataset and incorporate chi2 test in somatic test  
==============================

Purpose
------------

Different from previous somatic test where the percentage of posterior mu for case and control within an interval was computed and test, the new somatic test test on the mean of posterior mu for case and control sample. This helps because what we only need to set is just a threshold tau for case and  control independently rather than a interval and a threshold. 

In the new somatic test, combined the apriori knowledge of the tumor sample purity is 0.8, we set the threshold for homozygote mutation as 0.75. We set the threshold for reference (non-mutated) as 0.05)

We also incorporate the **chi2 test** for the test of non-reference mutation to improve the specificity.



Conclusions
-----------------
Combining somatic test and incorporating chi2 test for somatic test achieve very good result as can be seen in Table 1.

Background
----------------

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
RUn `MuBarPlot.py` to generate tables with correct setting in source code `rvd27.py`.

Results
-----------
The new somatic test without chi2 test is more sensitive than the previous somatic test and detected 15 likely-to-be-true somatic mutated positions as highlighted in blue in **Table 2**. However, there are 19 positions that are very likely to be false positives. The chi2 test will help removing the false positives as shown in **Table 1**. 

**Table 1 Positions called by new somatic test incorporating chi2 test**
![]('somaticcall_chi2.png')

**Table 2 Positions called by new somatic test without incorporating chi2 test**

![]('somaticcall_nochi2.png')

**Figure 1 The bar plot for mu in positions called by new somatic test incorporating chi2 test**
![]('HCC1187_mu_barplot.png')




Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
