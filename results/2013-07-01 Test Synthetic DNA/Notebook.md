2013-07-01 Test Synthetic DNA and Mutation Detection Based on RVD2.7 and Chi Square Test
==============================

Purpose
------------
To test RVD2.7 on real sequencing data where the truth about the mutation status is know and use chi square test to detect mutation position

Conclusions
-----------------
The chi square test failed to detect mutation with acceptable accuracy. 

Background
-----------------
A threshold of 1 was set to select positions with high potential to be a mutant position. For this positions, it was assumed that if the position was not a functional mutation, it should have even probability to mutate to the three non-reference bases, following multinomial distribution with probability vector [1/3, 1/3, 1/3]. To test whether these positions actually followed the assumed multinomial distribution, chi square test was performed for the goodness of fit. 

Materials and Equipment
------------------------------
Run `runall.py` first to get the hdf5 dataset, and then run `Test Synthetic DNA Notebook.py` to reproduce the analysis.

Experimental Protocol
---------------------------
Chi square test was consolidated in the previous RVD2.7 model. Parameter lambda have several options with default option as 2/3.

        lamda=1 Pearson's chi-square

        lamda=0 the log likelihood ratio statistic/ G^2

        lamda=-1/2 Freeman-Tukey's F^2

        lamda=-1  Neyman modified chi-square

        lamda=-2  modified G^2

About the goodness of fit test please go to [N. Cressie and T. R. C. Read, "Multinomial Goodness-of-Fit Tests," Journal of the Royal Statistical Society Series B-Methodological, vol. 46, pp. 440-464, 1984.](http://www.jstor.org/discover/10.2307/2345686?uid=3739696&uid=2&uid=4&uid=3739256&sid=21102472046131 "Goodness of fit") for reference


Results
-----------
For each dilution SNR and p value for chi square test were plotted. A eps=10^-30 was added to p value to avoid infinite values. For each positions there were six samples and chi square test was done for each position in each sample. For each position, log(pmax), log(pmin) and log(p) were plotted.
At the same time, sequencing result for positions with SNR>1 was provided for each dilution. 

![](http://i.imgur.com/S0X7EBC.jpg)![](http://i.imgur.com/r3gniT4.jpg)

![](http://i.imgur.com/dWB7oMe.jpg)![](http://i.imgur.com/6UgilLK.jpg)

![](http://i.imgur.com/HWcTaCN.jpg)![](http://i.imgur.com/YCd6rDC.jpg)

![](http://i.imgur.com/IXF2iHl.jpg)![](http://i.imgur.com/7hlM8M4.jpg)

![](http://i.imgur.com/k2hqd0e.jpg)![](http://i.imgur.com/IuUo1yM.jpg)

As can be seen from the figures, the chi square test didn't do well in mutation detection. For some unknown reason, some non-reference positions, the non-reference bases severely deviate from multinomial distribution with even probability for each category. Therefore, the assumption that the whole test depended on didn't hold water. 

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________________ Date: _____________________


Witnessed by: ________________________
