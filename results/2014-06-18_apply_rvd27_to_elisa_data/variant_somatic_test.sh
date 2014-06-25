#!/bin/sh
echo '----------------------'
echo 'Performing posterior density test'
RVD27=../../src/python/rvd27/rvd27.py
hdf5dir=./hdf5folder
controlsample=$hdf5dir/T0diploid_S2.hdf5
# casesample=$hdf5dir/T0haploid_S1.hdf5
# casesample=$hdf5dir/T1diploid_S3.hdf5
casesample=$hdf5dir/T2diploid_S4.hdf5
output=${casesample##*/}
output=${output%.*}

python $RVD27 somatic_test $controlsample $casesample -s 199096 -o $output