2013-11-25 Plot ROC with varying read depth by rvd29
==============================

Purpose
------------
To test the detection power of our graphical model using Bayesian Hypothesis Testing or including Chi-Square testing as well.

Conclusions
-----------------
For dilution rate as low as 0.1%, our graphical model using Bayesian Hypothesis Testing has high detection power when read depth is high enough. Including Chi-square Testing can limit False Positive Rate at a sacrifice of low True Positive Rate when read depth is low, but when read depth is high enough, the chi-square test does not decrease the True Positive Rate anymore.
### ROC shows the rvd29 with Jeffreys prior performs a little better than no prior situation especially on the small read depth.  

Background
-----------------
We use a two-stage test for variant at each position. First, we test whether the error rate in the case sample is significantly greater than the error rate in the control data for each position. Then, we test the hypothesis that the error rate is high is due to excess reads with a particular non-reference base. We accomplish the first test with a Bayesian posterior density test and the second with a chi-square goodness-of-fit test to a uniform distribution.

We tested the mutation detection power of our graphical model on synthetic dataset with dilution rates 0.1%, 0.3%, 1.0%, 10.0% and read depth varying from 10^2 to 10^5.Different read depth dataset was achieve by thining dataset with read depth up to 10^6 using Picard.

#### Bayesian Hypothesis Testing

To detect a mutation in a specific position, a threshold T is set for the posterior difference distribution testing.


 H1: muCase - muControl > T 
 
 H2: muCase - muControl <= T


We can calculate the posterior probability of the hypothesis H1 being true by integrating the posterior density over the correct region.

We accept the hypothesis H1 if this posterior probability is no lower than 95%, which means we are 95% sure that there is a mutation in this position.

ROC curve is plotted to evaluate the performance of the detection algorithm with different threshold for different mutation rate. 

#### Chi-Square Testing
We assume that at one position in one experimental replicate, the non-reference reads follows multinomial distribution of k = 3 categories, with probability vector p = (p1, p2, p3) where pi>0 for i =1, 2, 3 and p1+p2+p3=1.

The non-reference reads can be because of either mutation or various errors including sampling error, inter-replicated error. In theory, the non-reference read counts caused by non-oriented errors should be uniformly distributed among three non-reference categories, namely p_1=p_2=p_3=1/3. However, the distribution of non-reference read counts among three non-reference categories caused by mutation theoretically is biased. 

The null hypothesis can be denoted as

 H0: p=p0, where p0= (p1,p2,p3) and p1=p2=p3=1/3

The null hypothesis is tested using power-divergence family of statistics proposed by Cressie and Read (1984). Please refer to [Multinomial Goodness-of-Fit Tests](http://www.jstor.org/stable/2345686) for more information.

For each position we have six chi-square statistics from 6 experimental replicates. We use Fisher's method to combines the six statistics into a single one for inference. Bonferroni correction is applied to maintaining the familywise error rate (FWER). The significance level is chosen as 0.05.

It is optional whether to include Chi-square Testing for variant calling.If the read counts distribution across replicates accept the Bayesian hypothesis H1 and reject the chi-square hypothesis H0 in one position, a variant will be called in that position.


Materials and Equipment
------------------------------
run `subplotROC.py` to generate the ROC curve `ROC_subplots_with_chi2.eps`, set `chi2=True`

run `plotROCbyDilution.py` to generate the ROC curve `ROC_without_chi2.pdf`, set `chi2=False`

Experimental Protocol
---------------------------
Results
-----------

### Classification using Bayesian Hypothesis Testing

![](ROC_without_chi2_rvd29.png)
Figure 1.ROC curve varying read depth showing detection performance of model with Bayesian Hypothesis Test.

It can be seen that the model had high detection power when the read depth was at 10^5 level. At dilution rate at 0.1%, the model achieved performance with sensitivity at 99.74% and specificity at 92.86%. We missed one mutation at position 205 and had one false positive at position 8. The threshold was chosen at 0.022% which achieved the optimal sensitivity and specificity in the ROC curve. As expected, the detection power of the graphical model decreased when the read depth was reduced.

### Classification using Bayesian Hypothesis Testing and chi-square testing
![](ROC_subplots_with_chi2_rvd29.png)

Figure 2. ROC curve varying read depth showing method performance with Bayesian Hypothesis Testing and Chi-Square Test

Comparing to the detection results on using Bayesian Hypothesis Testing, there is no variant called or much fewer variants called when the read depth is low if including Chi-square testing. This means Chi-square testing is a stringent test when the read depth is low. Chi-square test limits **False Negative Rate** when the read depth is low (from 40 up to 4000), which also means we are confident in our variant callings if there is any. 

(However, if the read depth can be up to 40000, chi-square test still limits **False Negative Rate** though not strictly to zero, and not at a sacrifice of **True Positive Rate** anymore.)

One thing to be noted is that in chi-square test, if the sum of non-reference bases in one position is zero, the chi-squre statistics returned is NaN. NaN is set to fail to reject the null hypothesis.
Another thing is that to show the curves better, the xstick and ystick of the figures were limited to (-0.03,1.03)


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____ _Fan Zhang________        Date: ______ __11.25.2013_______


Witnessed by: ________________________
