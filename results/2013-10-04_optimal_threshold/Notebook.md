2013-10-04 optimize threshold for posterior difference bayesian hypothesis testing
==============================

Purpose
------------
To test the detection power of our graphical model using Bayesian Hypothesis Testing or including Chi-Square testing as well.

Conclusions
-----------------


Background
-----------------
we test whether the error rate in the case sample is significantly greater than the error rate in the control data for each position. We accomplish the first test with a Bayesian posterior density test.

We tested the mutation detection power of our graphical model on synthetic dataset with dilution rates 0.1%, 0.3%, 1.0%, 10.0% and read depth varying from 10^2 to 10^5.Different read depth dataset was achieve by thining dataset with read depth up to 10^6 using Picard.

#### Bayesian Hypothesis Testing

To detect a mutation in a specific position, a threshold T is set for the posterior difference distribution testing.


 H1: muCase - muControl > T 
 
 H2: muCase - muControl <= T


We can calculate the posterior probability of the hypothesis H1 being true by integrating the posterior density over the correct region.

We accept the hypothesis H1 if this posterior probability is no lower than 95%, which means we are 95% sure that there is a mutation in this position.






Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
Results
-----------
vcf files are available in current directory showing results with threshold optimized by Euclean distance, L1 distance and MCC respecitively.

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Yuting He________        Date: _________10/04/2013_________


Witnessed by: ________________________
