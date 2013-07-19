2013-07-17 Test Synthetic DNA T-test
==============================

Purpose
------------
Apply t-test on  data of mu obtained from MCMC sampling based on RVD2.7 model; real sequencing data where the truth about the mutation status is know was used.
empirical posterior distribution of

Conclusions
-----------------
This experiment showed certain detection ability of the approach based on rvd2.7 and follow-up t-test, pretty much similar to the results achieved with SNR method. The rationale of this method is not strictly sound.  

Background
-----------------
This experiment is about applying unpaired  t test to compare the means of control mu samples and case mu samples for each position. About independent two-sample t-test, please go to [Independent two sample t-test](https://en.wikipedia.org/wiki/Student's_t-test#Independent_two-sample_t-test) for reference.

There are several things to be considered in this experiment.

1) The prior of mu is beta distribution, while the posterior distribution of mu is hard to find out
   T-test assumes that the data is from a normal distribution; however, it is unknown to all what's the posterior distribution of mu is. 
   Normal approximation to the Beta distribution is proper when (alpha+1)/(alpha-1).=1, and (beta+1)/(beta-1).=1. A rule of thumb is that alpha and beta are both equal to 10 or more. 

2) Assuming the posterior distribution of mu is or can be approximated as normal distribution so independent two sample t-test can be done, is it reasonable to assume the variance of mu for control and case data for one specific position are equal? 

3) t test is based on frequentist statistics. Performing t-test on empirical posterior distribution of mu is improper. 

Experiment was done and results are as shown below.

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------


Results
-----------
Figure 1(a) ![](http://i.imgur.com/q9Lmh9g.jpg)

Figure 1(b) ![](http://i.imgur.com/EyfaO0L.jpg)

Figure 1(c) ![](http://i.imgur.com/93mKPQR.jpg)

Figure 1(d) ![](http://i.imgur.com/PhDF65p.jpg)

Figure 1(e) ![](http://i.imgur.com/VmbcXwQ.jpg)

Figure 1 are histograms of sampling data of mu for two positions, non-mutated position 50 and mutated position 245. As can be see from the figures, 1) the sampling data of mu are roughly normally distributed; 2) for all dilution, the case and control mu are much more distantly distributed in mutated position 245 than in non-mutated position 50, which suggests the feasibility of t-test. 

![](http://i.imgur.com/jdDSWOj.jpg) 

Figure 2 t value when assuming the variance of case mu and control mu are equal

![](http://i.imgur.com/HWy8Fkw.jpg)

Figure 3 t value when assuming the variance of case mu and control mu are unequal

From Figure 2 and Figure 3 it can be seen that the t-test approach achieved similar results as the SNR approach. Under each dilution, the t value of mutated positions are generally distinctively higher than non-mutated positions. However, there are mutated positions have lower t value than non-mutated positions, which affects the overall mutation detection accuracy. 

Comparing Figure 2 and Figure 3, the difference between results assuming equal variance and unequal variance are not obvious. 


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ______Yuting He__________ Date: _____________________


Witnessed by: ________________________
