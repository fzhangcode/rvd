2013-11-25 Read depth and precision parameter M evaluation across position
==============================

Purpose
------------
Evaluate the of properties of precision parameter M across positions. M is using Jeffreys' prior.

Conclusions
-----------------
The precision parameter is low in the two ends of positions, while is fairly high in the middle positions. 
###The M varies across positions. A important observation from dilution 100.0% and 10.0% is that the M value of the mutation locations are smaller than 100, which isn't seen in Gammar prior situation. For non-mutation locations M are generally in between 10^4 and 10^5.
### Additionally, the distribution of M shows more normal and stable than the M without priors.

Background
-----------------


Materials and Equipment
------------------------------
Data: hdf5 data sets in folder 2013-11-18_ test_ rvd29_ synthetic_ data_ dc=10

Codes: M_loc.py 

Experimental Protocol
---------------------------


Results
-----------
![](M_loc.pdf)


Figure 1. Precision parameter Mj across positions. The y axis in left panels are in linear scale, while the y axis in the right panels are in log scale.

From the figure it can been seen that the precision parameter is low in the two ends of positions, while is fairly high in the middle positions. The M varies across positions, but are generally in between 10^4 and 10^5 for most positions. A important observation from dilution 100.0% and 10.0% is that the M value for the mutation location are smaller than 100.

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: ________Fan  Zhang_______     Date: __________11.25.2013 ___________


Witnessed by: ________________________
