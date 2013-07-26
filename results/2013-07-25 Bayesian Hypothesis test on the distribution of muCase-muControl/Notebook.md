2013-07-20 Test on Metropolis Hasting sampling size with Synthetic DNA
==============================

Purpose
------------
To test the effect of Metropolis Hasting sampling sampling size on mutation detection

Conclusions
-----------------
The Metropolis Hasting sampling sampling size has noticeable effect on t value. The optimal value is not yet found since the program is till running.


Background
-----------------
Given that Gibbs sampling size 400 is enough to build the empirical posterior distribution of mu_s, Metropolis Hasting sampling size of 1,10,50,100,1000 were tried to test the effect of Metropolis Hasting sampling size on mutation detection. 

Materials and Equipment
------------------------------


Experimental Protocol
---------------------------


Results
-----------
![](http://i.imgur.com/rJYDLAD.jpg)
Figure 1 t-test statistic when gibbs sampling size is 400 and MH sampling size is 1
![](http://i.imgur.com/fDPKzay.jpg)
Figure 2 t-test statistic when gibbs sampling size is 400 and MH sampling size is 10
![](http://i.imgur.com/NcjXDLH.jpg)
Figure 3 t-test statistic when gibbs sampling size is 400 and MH sampling size is 50
![](http://i.imgur.com/CrsrVu1.jpg)
Figure 4 t-test statistic when gibbs sampling size is 4000 and MH sampling size is 50
![](http://i.imgur.com/viyIXGz.jpg)
Figure 5 t-test statistic when gibbs sampling size is 400 and MH sampling size is 100
![](http://i.imgur.com/f2qYVxA.jpg)
Figure 6 t-test statistic when gibbs sampling size is 400 and MH sampling size is 1000

Several things to mention:
1. Figure 4 is from the test on 2013-07-18. It's included here to compare with Figure 3. In the two figures, the MH sampling size is the same,50, while the gibbs step size is different, 400 and 4000. From the figure it can be seen that the difference is very insignificant, which further prove that gibbs sampling size at 400 is enough to build a full posterior of mu.

2. Comparing all the figures, it can be observed, though not very obvious and not quantified yet, that the MH sampling size at 1000 separate the mutation positions and the non-mutation positions best. 

3. Till now it can be seen the metropolis hasting sampling size 100 achieved best performance. To find out the optimal/sufficient sampling size, program with sampling size of 1000,5000 is under running. At dilution 0.1, except position 205, all other mutated positions have higher t-value than non-mutated positions. at dilution 0.3, the 4 mutated positions are not extinct from some extreme non-mutated positions. Result at 0.3 is the worst among all five dilutions.For dilution at 1.0 and dilution 10.0, the mutated position are all distinctive from others. While for dilution at 100, there are some high t test statistic values in position at really close 400. This phenomenon happened because the last part of gene has really low read depth compare to other positions. 

4. It is noticed that in all figures, the result in dilution at .1 is better than that in dilution at .3, why?

5. we noticed that there is a very strange point at position 205 in dilution 0.1. By looking in to the sequence data as shown in Table 1, it can be seen that the distribution of the data in position 205 is not as extreme as mutated position 145, but similarly as extreme as none position 146. We need to think about whether this is because of the data or the insensitive of the model that can not distinct position 205 from non-mutated positions.

Table 1

![](http://i.imgur.com/ceosv9o.jpg)



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: __________Yuting He______  Date: _____________________


Witnessed by: ________________________
