2013-07-25 Bayesian Hypothesis test on Z posterior distribution 
==============================

Purpose
------------
Bayesian Hypothesis test and ROC analysis on Z posterior distribution.

Conclusions
-----------------
The ROC plots were very good, indicating the Bayesian Hypothesis testing can be very promising at rare mutation detection.

A issue which needs to solve is how to find the optimal threshold. 

Background
-----------------
This experiment aimed at comparing the difference between muControl_s and muCase_s. Define X as the simulated posterior distribution of muCase_s, and Y as the simulated posterior distribution of muControl_s, the difference posterior distribution is defined as Z=X-Y. The Z posterior distribution was obtained by doing N times of following procedures: randomly get x from X; randomly get y from Y;z,one sample of Z, is z=x-y.

With Z distribution obtained, one sided Bayesian Hypothesis was done by setting an threshold TH. with significance alpha=0.05, p(Z>TH) was computed and compared with 1-alpha to classify the position.

ROC analysis was done in the experiment to evaluate the performance of the classification method. ROC analysis also gave thoughts on how to determine threshold.

In this experiment, the gibbs sampling size is 4000, MH sampling size is 50.The data used in this experiment was synthetic dataset with known mutation positions.


Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
run runall.py to get the hdf5 dataset

run BayesianHypoTest.py to reproduce the result for Bayesian Hypothesis Test and ROC analysis


Results
-----------
Issue to resolve: How to determine the Threshold?
![ROC plot](ROCplot_nGibbs=4000_nMH=50_steps=100_alpha=0.05.png)

Dilution=0.1
![](histXY_nGibbs=4000_nMH=50_dilution=0_1.png)
![](histZ_nGibbs=4000_nMH=50_dilution=0_1.png)


Dilution=0.3
![](histXY_nGibbs=4000_nMH=50_dilution=0_3.png)
![](histZ_nGibbs=4000_nMH=50_dilution=0_3.png)


Dilution=1
![](histXY_nGibbs=4000_nMH=50_dilution=1_0.png)
![](histZ_nGibbs=4000_nMH=50_dilution=1_0.png)

Dilution=10
![](histXY_nGibbs=4000_nMH=50_dilution=10_0.png)
![](histZ_nGibbs=4000_nMH=50_dilution=10_0.png)

Dilution=100
![](histXY_nGibbs=4000_nMH=50_dilution=100_0.png)
![](histZ_nGibbs=4000_nMH=50_dilution=100_0.png)



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: __________Yuting He______  Date: _____________________


Witnessed by: ________________________
