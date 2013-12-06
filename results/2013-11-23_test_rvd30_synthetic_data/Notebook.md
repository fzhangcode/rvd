2013-11-25 Test rvd30 on synthetic data
==============================

Purpose
------------
1. Test rvd30 the lognormal prior performance
2. Read depth and precision parameter M evaluation across position
3. To test the detection power of our graphical model using Bayesian Hypothesis Testing or including Chi-Square testing as well.

Conclusions
-----------------
From figure2, we know lognormal prior performs very good especially on the low read depth, better than the tests before.

Background
-----------------
Refer to rvd29 folder background.


Materials and Equipment
------------------------------
1.  1) run "runall.py" to acquire hdf5 files, 
 
    2) run "test_rvd30" to test the results, 
 
    3) run "character.py" to generate the results table.

2. run "M_loc.py" to plot the M.

3. run "plotROCbyDilution.py" to generate the ROC figure.

Experimental Protocol
---------------------------
Results
-----------
###1. See table: "character_rvd30.xls"

###2. Read depth and precision parameter M evaluation across position

![](M_loc_rvd30.png)
Figure 1. Precision parameter Mj across positions. The y axis in left panels are in linear scale, while the y axis in the right panels are in log scale.

From the figure it can been seen that the precision parameter is low in the two ends of positions, while is fairly high in the middle positions. The M varies across positions, but are generally in between 10^4 and 10^5 for most positions. A important observation from dilution 100.0% and 10.0% is that the M value for the mutation location are smaller than 100.

###3. Classification using Bayesian Hypothesis Testing

![](ROC_without_chi2_rvd30.png)
Figure 2. ROC curve varying read depth showing detection performance of model with Bayesian Hypothesis Test.

It can be seen that the model had high detection power when the read depth was at 10^5 level. At dilution rate at 0.1%, the model achieved performance with sensitivity at XXX and specificity at XXX. We missed one mutation at position 205(?) and had one false positive at position 8(?). As expected, the detection power of the graphical model decreased when the read depth was reduced.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____ _Fan Zhang________        Date: _____ __11.25.2013_______


Witnessed by: ________________________
