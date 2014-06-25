#!/bin/sh
echo '----------------------'
echo 'Performing posterior density test'
RVD27=../../src/python/rvd27/rvd27.py
hdf5dir=./hdf5folder
controlsample=$hdf5dir/T0diploid_S2.hdf5

output=${controlsample##*/}
output=${output%.*}

python $RVD27 germline_test $controlsample -o $output -i 0.01 1