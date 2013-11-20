2013-09-21 SNP calling using RVD2 with threshold set at zero
==============================

Purpose
------------
To evaluate the variant detection power of RVD2 with threshold set at zero

Conclusions
-----------------
RVD2 with threshold set at zeros shows high variant detection power.

Background
----------------

Materials and Equipment
------------------------------
halfdilution.py


Experimental Protocol
---------------------------

Results
-----------

The program generates vcf files which saves variant calling results. They are available in the same directory.

Using the characer.py in `../2013-09-19_operating_characteristics` it is possible to generate a summarizing table showing the sensitivity/specificity for variant calling using RVD2 with threshold set at zero. A screenshot is shown as follows:

![](http://i.imgur.com/VtQIIEZ.png)

Figure 1. Sensitivity/specificity achieved using  RVD2 with threshold set at zero for variant calling.

It can be seen that RVD2 with threshold set at zero shows high variant detection power across many dilution and coverage combination. rvd2_zero achieved sensitivity/specificity=1.00/1.00 for all read depth when dilution rate at 100% and rvd2_zero achieved sensitity/specificity=1.00/1.00 for dilution rate at 10% except when read depth is too low at 22. Coverage at 535 is enough for rvd2_zero to achieve sensitivity/specificity=0.93/0.95 when dilution rate is 1.0%. For dilution rate at 0.1%, sensitity/specificity=0.86/0.98 when read depth is as high as 41449.


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____Yuting He____________     Date: _______09/30/13______________


Witnessed by: ________________________
