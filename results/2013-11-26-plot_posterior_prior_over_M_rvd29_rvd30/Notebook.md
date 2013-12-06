2013-12-02 Plot the distribution of priors and posteriors for rvd29 (Jeffreys prior) and rvd30 (lognormal prior)
==============================

Purpose
------------
Evaluate the performance of priors and posteriors of Mj over different M values. We want to compare the two priors: Jeffreys prior from rvd29, and lognormal prior from rvd30.

Conclusions
-----------------
From the figures, The posterior distribution shows normal and stable. It doesn't show a peak at some unreasonable point or a strange shape. 

Background
-----------------
To investigate the performance of different priors, the assess from the posterior view should be taken. We are not going to compare which prior is better by judging posterior, but we can analysis the results, and research the sensitivity of the model to the choice of the prior from the shape of the posterior distribution. 


Materials and Equipment
------------------------------
To plot the Jeffreys and lognormal priors, the code is the same as the relevant code in rvd29 and rvd30. 

To plot the Jeffreys and lognormal posterior, the Gaussian Kernel Density Estimate is applied.

Finally we normalize the figures to make them shown in the same scale.



Data: The datasets are under the same experiment condition: read depth rate =1/100, and to acquire more samples here "gibbs_sample=20000" in gibbs sampling, instead of 400 samples.
  
    Jeffreys prior hdf5 data sets are in folder ./hdf5_rvd29_gibbs=20000_dc=100

    Lognormal prior hdf5 data sets are in folder ./hdf5_rvd30_gibbs=20000_dc=100

Codes: plot\_prior\_posterior.py 

(In addition, "plot\_prior\_posterior\_rvd29.py" and "plot\_prior\_posterior\_rvd30.py" are respectively the codes for plotting rvd29 and rvd30.)

Experimental Protocol
---------------------------
Mutant and non-mutant locations are chosen to see the probability distribution for the priors and posteriors. The chosen locations are from the around the middle place of the base length. 

Results
-----------
The figure covers different locations (mutation and non-mutation) to compare the performances. The mutations are chosen in the middle of the dase length.


![](prior_posterior_dilu=10.png)
Figure 1. Probability distribution of priors and posteriors, when dilution=10%

![](prior_posterior_dilu=0.3.png)
Figure 2. Probability distribution of priors and posteriors, when dilution=0.3%

![](prior_posterior_dilu=1.0.png)
Figure 3. Probability distribution of priors and posteriors, when dilution=1.0%

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Fan  Zhang__ _____     Date: ________12.03.2013 _________


Witnessed by: ________________________
