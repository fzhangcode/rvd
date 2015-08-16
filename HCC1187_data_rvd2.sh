#!/bin/sh
echo '----------------------'
echo 'Running RVD2 in the HCC1187 subject PAXIP1 gene dataset.'
echo 'Please provide the directory to the control and case hdf5 files.'
echo '----------------------'
echo 'Generating control and case hdf5 file from depth chart files:'
python rvd27.py gibbs $1 -o control_HCC1187 -p 10 -g 4000 -m 5 -s 199096
python rvd27.py gibbs $2 -o case_HCC1187 -p 10 -g 4000 -m 5
# echo '----------------------'
echo 'Performing hypothesis test (germline test and somatic test)'
echo '----------------------'
python rvd27.py germline_test control_HCC1187.hdf5 
echo '----------------------'
python rvd27.py somatic_test control_HCC1187.hdf5 case_HCC1187.hdf5
echo 'Done.'