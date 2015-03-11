#!/bin/sh
#vpath %.dc ./depth_chart/3000/
vpath %.dc ./depth_chart/10/

RVD2=../../src/python/rvd27/rvd27.py

fit: 20100916_c1_p2.04_ATT_1000.hdf5 20100916_c1_p2.04_ATT_2000.hdf5 20100916_c1_p2.04_ATT_3000.hdf5 20100916_c1_p2.04_ATT_4000.hdf5 

# Fit the RVD2 model

#20100916_c1_p2.04_ATT_?.hdf5 : 20100916_c1_p2.04_ATT.dc
#	time python $(RVD2) gibbs $^ -o $(basename $@) -p 1 -g 1000 -m 5 -s 19860522


20100916_c1_p2.04_ATT_1000.hdf5 : 20100916_c1_p2.04_ATT_1000.dc
	time python $(RVD2) gibbs $^ -o $(basename $@) -p 1 -g 1000 -m 5 -s 19860522 > time_1000.txt

20100916_c1_p2.04_ATT_2000.hdf5 : 20100916_c1_p2.04_ATT_2000.dc
	time python $(RVD2) gibbs $^ -o $(basename $@) -p 1 -g 2000 -m 5 -s 19860522 > time_2000.txt

20100916_c1_p2.04_ATT_3000.hdf5 : 20100916_c1_p2.04_ATT_3000.dc
	time python $(RVD2) gibbs $^ -o $(basename $@) -p 1 -g 3000 -m 5 -s 19860522 > time_3000.txt

20100916_c1_p2.04_ATT_4000.hdf5 : 20100916_c1_p2.04_ATT_4000.dc
	time python $(RVD2) gibbs $^ -o $(basename $@) -p 1 -g 4000 -m 5 -s 19860522 > time_4000.txt

