#!/bin/sh
echo '----------------------'
echo 'Generating hdf5 files'
RVD27=../../src/python/rvd27/rvd27.py

dcfile=./depth_chart/T0diploid_S2.dc

mkdir -p hdf5folder
hdf5file=${dcfile##*/}
hdf5file=${hdf5file%.*}
sufix='_newseed'
hdf5file=$hdf5file$sufix
hdf5file=./hdf5folder/$hdf5file
echo $hdf5file

python $RVD27 gibbs $dcfile -o $hdf5file -p 10 -g 4000 -m 5 -s 19891129