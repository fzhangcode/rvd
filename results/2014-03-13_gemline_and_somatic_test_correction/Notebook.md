2014-03-13 Correct Germline and Somatic mutation detection method 
==============================

Purpose
------------
Corrected Germline and Somatic mutation detection method echos to the result section in the paper.

Conclusions
-----------------

Background
----------------

We used RVD2 to identify germline and somatic mutations in the diploid HCC1187 sample. To identify germline mutations, we compute the empirical mean of posterior control mu. If the posterior mean falls between [0, 0.05), the position is considered as homozygous reference; if the posterior mean falls between [0.5, 0.0.75), the position is considered as heterozygous mutant; otherwise, if the posterior mean falls between  [0.75, 1.0], the position is homozygous mutation.

To identify somatic mutations, we considered scenarios when the case(tumor) error rate is lower than the control(germline) error rate (e.g. loss-of-heterozygosity) as well as scenarios when the case(tumor) error rate is higher than the control(germline) error rate (e.g. homozygous somatic mutation). The two hypothesis tests are then $\Pr( \mu_j^{\Delta} \geq \tau ) > 1-\alpha$ and $\Pr( \mu_j^{\Delta} \leq \tau ) > 1-\alpha$. The size of the test is $\alpha=0.05$.




Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
RUn `MuBarPlot.py` to generate tables with correct setting in source code `rvd27.py`.

Results
-----------
![]('HCC1187_mu.png')

**Figure 1 The bar plot for mu in positions called by new somatic test and Germline test. Positions with star (\*) were called by Somatic test, where case mu and control mu are significantly different from each other.**



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 02/24/14


Witnessed by: ________________________
