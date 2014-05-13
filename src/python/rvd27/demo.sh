#!/bin/sh
echo '----------------------'
echo 'RVD2 demo.'
echo '----------------------'
echo 'Generating simulation data:'
python rvd27.py gen -s 199096
echo '----------------------'
echo 'Performing hypothesis test (paired_difference_test)'
echo '----------------------'
python rvd27.py paired_difference_test control_simulation.hdf5 case_simulation.hdf5 
echo 'Done. Succeed if a variants_paired_difference.hdf5 file and a variants_paired_difference.vcf file are created.'