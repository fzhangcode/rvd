2013-10-04 surface fitting for optimal threshold
==============================

Purpose
------------
do surface fitting to find the optimal threshold relationship with dilution rate and coverage

Conclusions
-----------------
run 

Background
-----------------
run SNP_fitT.sh to generate final vcf files.

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------
Second order polynomial surface fitting was done to optimal threshold generated according to MCC.
 
Results
-----------

**Surface fitting**

x represents dilution rate(0.1,0.01,0.001,0.0001, while 1 is excluded); y represents median coverage(41342, 4147, 403, 37)


 Linear model Poly22:  Badly conditioned.
     fitresult1(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
     Coefficients (with 95% confidence bounds):
       p00 =   -0.004613  (-0.008548, -0.0006773)
       p10 =      0.2574  (-0.3409, 0.8556)
       p01 =   5.926e-07  (-6.389e-07, 1.824e-06)
       p20 =      -1.325  (-7.062, 4.411)
       p11 =  -2.323e-06  (-4.801e-06, 1.549e-07)
       p02 =  -1.114e-11  (-4e-11, 1.772e-11)



 Linear model Poly22:

x represents log10(dilution rate), y represents log10(coverage)
     fitresult2(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2
     Coefficients (with 95% confidence bounds):
       p00 =      0.0246  (0.007546, 0.04166)
       p10 =     0.02307  (0.01126, 0.03488)
       p01 =    0.001648  (-0.00631, 0.009606)
       p20 =    0.002653  (-9.827e-05, 0.005405)
       p11 =   -0.002472  (-0.003908, -0.001036)
       p02 =  -0.0008988  (-0.002074, 0.0002759)

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Yuting He________        Date: _________10/04/2013_________


Witnessed by: ________________________
