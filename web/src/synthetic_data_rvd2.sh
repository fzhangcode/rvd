#!/bin/sh
set -e 
echo '----------------------'
echo 'Running RVD2 in the synthetic dataset.'
echo '----------------------'
echo 'Please enter the full  directory to depth chart files for control sample:'
read control_dir
python rvd27.py gibbs $control_dir -o control_synthetic -p 10 -g 4000 -m 5 -s 199096
echo 'Please enter the full  directory to depth chart files for case sample:'
read case_dir
python rvd27.py gibbs $case_dir -o case_synthetic -p 10 -g 4000 -m 5
echo '----------------------'
echo 'Performing hypothesis test (germline test and somatic test)'
echo '----------------------'
python rvd27.py paired_difference_test control_synthetic.hdf5 case_synthetic.hdf5
echo 'Done.'