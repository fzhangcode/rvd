2015-03-10 timing for rvd2 with different ngibbs 
==============================

Purpose
---
Measure the timing for rvd2 model.


Conclusions
-----------------
The time increases when the depth chart increases. But the time is independent with the median depth chart.

Background
----------------



Materials and Equipment
------------------------------
Pool in gibbs function is used specifically in sampleMuMh function.  
Sample Mu across different locations using parallel.
Region size: 400 bps(1-400)

dc=130:  
3000/20100916_c1_p2.04_ATT.dc is from 2015-02-23_run_rvd2_SS_FDR_downsample_depth_100\depth_chart\3000\ 

dc=40000  
10/20100916_c1_p2.04_ATT.dc is from Y:\yhe2\Research\rvd2\results\2013-08-06_Downsample_Read_Depth\depth_chart\10

Experimental Protocol
---------------------------
`make -j4`

Results
-----------
A timing table is made based on these two results.

For dc=130, 4 jobs run parallel:

	real    13m34.826s
	user    14m40.172s
	sys     0m3.316s
	
	real    26m30.395s
	user    28m36.953s
	sys     0m11.873s
	
	real    39m54.340s
	user    42m53.443s
	sys     0m29.447s
	
	real    55m6.064s
	user    58m48.899s
	sys     0m41.115s

For dc=40000, 4 jobs run parallel:
	real    13m0.044s
	user    14m5.003s
	sys     0m3.860s
	
	real    26m13.154s
	user    28m19.301s
	sys     0m11.940s
	
	real    40m22.544s
	user    43m27.816s
	sys     0m23.675s
	
	real    53m42.820s
	user    57m32.277s
	sys     0m40.367s

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Fan Zhang     Date: 2015-03-10


Witnessed by: ________________________
