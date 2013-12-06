2013-11-12 rvd 28: Gamma prior on Mj
==============================

Purpose
------------
Build Gamma prior for Mj.

Conclusions
-----------------
Maybe Gamma prior cannot give a good estimate for shape parameter a.

Background
-----------------
The initial values for the model parameters and latent variables is obtained
by a method-of-moments (MoM) procedure. MoM works by setting the population moment equal to the sample moment.Sampling the posterior over Mj. (see 2013-11-15_Priors_Manuscript.pdf for detail inference process)

Materials and Equipment
------------------------------
Apply Gamma prior to Mj on the model of rvd27.

Experimental Protocol
---------------------------
Run rvd28\_test.py, use generate_sample and estimate_mom to test the rvd28 model.

Results
-----------
Compare the parameters a, b, mu0 and M0 before generate_sample and after estimate_mom, the input we choose are: a=2.5e4, b=5, mu0=0.25, M0=10; the b, mu0 and M0 are very close by estimate_mom. But the parameter a is far away. 


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ____  Fan Zhang________     Date: ____2013/11/15 __


Witnessed by: ________________________
