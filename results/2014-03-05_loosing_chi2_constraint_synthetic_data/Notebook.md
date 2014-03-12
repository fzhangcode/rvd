2014-03-04 loosing chi2 p-value constraint on synthetic data
==============================

Purpose
------------
As the chi2 test (with Bonferroni correction) is too strict for HCC1187 test, this experiment aims at testing how loosing chi2 p-value constraint will effect ss/fdr in synthetic dataset.

Conclusions
-----------------
Combing the Sensitivity and Specificity table and False Discovery table, the "chi2, mu0=1" option, which means  chi2 p-value is 0.05 without Bonferroni correction and subtracting mu0 in the bayesian hypothesis test achieve the most preferable result.

Background
----------------

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------

Results
-----------
Please see the histogram for this positions. 
![]('SS.png')
From the view of Sensitivity and Specificity, 
![]('FDR.png')

From the table it can be seen that loosing chi2 p-value will improve sensitivity at sacrifice of specificity. Also, the FDR will be degraded when chi2 p-value is loosed. 

Combing the Sensitivity and Specificity table and False Discovery table, the "chi2/J, mu0=1" option, which means  chi2 p-value is 0.05 with Bonferroni correction and subtracting mu0 in the bayesian hypothesis test achieve the most preferable result.

Regarding to the clinical dataset HCC1187, considering the long sequence length (44K) and low read depth(~50), it will be more preferable if we don't do Bonferroni correction. 

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
