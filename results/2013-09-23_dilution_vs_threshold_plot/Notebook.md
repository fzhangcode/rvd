2013-09-23 threshold vs dilution and coverage plot
==============================

Purpose
------------
To find the relationship between optimal threshold and dilution rate/coverage

Conclusions
-----------------
The three dimensional calibration plot for optimal threshold vs dilution rate and average coverage is generated in Figure 3.

Background
----------------

Materials and Equipment
------------------------------
hdf5 data files from folder :

	'2013-08-14_Compute_ROC_Synthetic_avg10',
	'2013-08-14_Compute_ROC_Synthetic_avg100',
	'2013-08-14_Compute_ROC_Synthetic_avg1000',
	'2013-08-14_Compute_ROC_Synthetic_avg10000'.

Codes:
`plotROCbyDilution_linear.py` and `plotROCbyDilution_log.py` under the same directory.

Experimental Protocol
---------------------------
1. Import Gibbs sampling data from the dataset file folders. 

1. Sample the difference posterior distribution 15 times (post_diff=Mucase-Mucontrol) with sampling size at N=1000. 

1. Generate the ROC curves and get the optimal threshold. 

1. Plot the optimal threshold across dilution rate and coverage. Spline is applied for in the plot optimal threshold vs dilution rate and average coverage. 


Results
-----------
![](OptimalT_vs_Dilution4_log.png)

Figure 1. Calibration plot for optimal threshold across different dilution rate (0.1%, 0.3%, 1.0%, 10.0%) with subplots differing in average coverage.

In each subplot in Figure 1, dilution rate axis is in log scale, and the optimal threshold axis is in linear scale. Average across 15 times of sampling is denoted in the plot, and +/-1 standard deviation is denoted by the red vertical bar in the plot. It can be seen expect when the dilution rate and coverage are both low ==(dilution rate=0.1 or 0.3, coverage average=42.0), the standard deviation is relatively high, otherwise, the standard deviation is insignificant comparing to the average. Therefore, the sample size N=1000 is large enough to get the empirical difference posterior distribution and the optimal threshold.


![](OptimalT_vs_Dilution5_log.png)

Figure 2. Calibration plot for optimal threshold across different dilution rate (0.1%, 0.3%, 1.0%, 10.0%, 100.0%) with subplots differing in average coverage.

The difference between Figure 1 and Figure 2 is that Figure 2 includes dilution rate at 100%. It can been seen that for the first four dilution rate, the optimal threshold are generally increasing. However, there is a sharp drop of optimal threshold from 10% to 100% dilution rate. This is because at the dilution 100%, there is a wide range of choice for threshold that can achieve optimal result (sensitivity=1, specificity=1), and the program returns the lowest threshold that can meet the requirement. That is to say, at some point between 10% and 100% dilution rate, we have wide choice of dilution rate. For convenience, we can set the optimal threshold at zero when the dilution rate is in range of 10.0% up to 100.0%.

![](OptimalT_vs_Dilution4_log_3D.png)
Figure 3. Three dimensional Splined calibration plot for optimal threshold across different dilution rate (0.1%, 0.3%, 1.0%, 10.0%) with subplots differing in average coverage.
 
The dilution rate axis and the coverage axis are in log10 scale.

It can be seen from Figure3 that in log scale,  optimal threshold is relative more stable across dilution rate than coverage. When the coverage is low, the optimal threshold has a decreasing trend in the begining, and then has a increasing trend. The ununiform trend is mainly introduce by 400, the optimal threshold stably increases as the dilution rate increase. Spline is applied in the three dimensional plot.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Yuting He________     Date: __________09/29/13__________


Witnessed by: ________________________
