2013-11-12 rvd 29: Jeffreys prior on Mj
==============================

Purpose
------------
Build a uniformative prior for Mj.

Conclusions
-----------------
Jeffreys priors give the same estimate for mu0 and M0 with Method of Moments.

Background
-----------------
Though conjugate priors are computationally nice, objective Bayesians instead prefer priors which do not strongly influence the posterior distribution. Such a prior is called an uniformative prior. 
Jeffreys priors are a generalization of these ideas, and can deliver a broad range of priors that incorporates these special cases. They are quite resonable in one dimension. They are based on a principle of invariance: one should be able to apply these priors to certain situations, apply a change of variable, and still get the same answer. [Michael I. Jordan's lecture]

Materials and Equipment
------------------------------
Apply Jeffreys prior to Mj on the model of rvd27.

Experimental Protocol
---------------------------
Run rvd29\_test.py, use generate_sample and estimate_mom to test the rvd29 model.

Results
-----------
Compare mu0 and M0 before generate_sample and after estimate_mom, the input we choose are: mu0=0.25, M0=10; the output by estimate_mom are around: 0.25376666, and 7.6696 by several tests. They are very close.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ____  Fan Zhang________     Date: ______2013/11/15  __


Witnessed by: ________________________
