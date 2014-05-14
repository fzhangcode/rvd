#!/bin/sh
set -e 
echo '----------------------'
echo 'Running RVD2 in the HCC1187 subject PAXIP1 gene dataset.'
echo '----------------------'
echo 'Please enter the full  directory to depth chart files for control sample:'
read control_dir
python rvd27.py gibbs $control_dir -o control_HCC1187 -p 10 -g 4000 -m 5 -s 199096
echo 'Please enter the full  directory to depth chart files for case sample:'
read case_dir
python rvd27.py gibbs $case_dir -o case_HCC1187 -p 10 -g 4000 -m 5
# echo '----------------------'
echo 'Performing hypothesis test (germline test and somatic test)'
echo '----------------------'
python rvd27.py germline_test control_HCC1187.hdf5 
echo '----------------------'
python rvd27.py somatic_test control_HCC1187.hdf5 case_HCC1187.hdf5
echo 'Done.'