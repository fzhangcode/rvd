2013-08-13 Compare ROC by Read Depth
==============================

Purpose
------------
To measure sensitivity and specificity for full RVD2 pipeline for various read depths on synthetic gene data.

Conclusions
-----------------
The chi-square test lowered the quality of ROC curves. 

Background
-----------------
Materials and Equipment
------------------------------
Experimental Protocol
---------------------------
Results
-----------

### result with chi-square testing
In chi-square test, the if the sum of non-reference bases in one position is zero, the chi-squre statistics returned is NaN. For one position if there is a NAN, the testing result will be purely based on Bayesian Hypothesis testing as well. To put it in programming language, I adapted the codes to:

`if postP[i] >0.95 and (chi2P[i] < 0.05/J or np.isnan(chi2P[i]))`
 
 `...`
 
![](ROC_with_chi2.png)

Figure 1. ROC curves from Bayesian Hypothesis Testing with chi-square test 


### result without chi-square testing

![](ROC_without_chi2.png)
Figure 2.ROC curves from Bayesian Hypothesis Testing without chi-square test


Up to now I only have generated two dataset, as used in this experiment. 
Comparing two figures it can be seen

- Chi2 test limited **False Negative Rate**, since even when all positions passed Bayesian Hypothesis Testing, some of them will fail Chi2 test. So that why the ROC curves in Figure 1 can not go up to point (1,1)
- The Chi2 test didn't improve **True Positive Rate**. On the contrary, in dilution 0.1 and 0.3, the optimal**True Positive Rate** is severely affected by the chi-square test. 
- To show the curves better, the xstick and ystick of the figures were limited to (-0.01,1.01)


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________     Date: _____________________


Witnessed by: ________________________
