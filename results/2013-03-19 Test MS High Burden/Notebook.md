2013-03-29 MS High Burden
==============================

Purpose
------------
To test RVD2.6 on actual clinical sequence data.

Conclusions
-----------------
RVD2.6 does not handle widely varying error rates well due to the M parameters tied across positions. I modified the model and updated the version to 2.7 to put an improper prior on M. Further development on a good prior may be needed.

Background
-----------------

Materials and Equipment
------------------------------

Experimental Protocol
---------------------------

Results
-----------
When running RVD2.6 on the data from CD226, I observed that the baseline error  rate would move to 0.30 - much higher than it should be. By eye, the r/n ratio for most positions is near 0.1%. Upon close inspection of a few positions, I found that while the vast majority of the error rates are 0.1%, some positions are heterozygous. So the "error rate" approaches 50%. That causes a huge difference compared to the mu0 expected. The algorithm adjusts the a and b to be extremely diffuse which in turn allows mu to match the overall average and essentially not contain any position specific information. 

I fixed M at the MoM estimate and found that the model worked very well on this test data set. My guess is that by tying these parameters across positions through the gamma prior, they were not able to float to match such wide differences between error rates.

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _________________________ Date: _____________________


Witnessed by: ________________________
