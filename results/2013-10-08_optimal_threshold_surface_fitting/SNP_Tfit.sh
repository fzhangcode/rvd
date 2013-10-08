#!/bin/sh

echo "Generate optimal threshold according to Euclidean distance, L1 distance and MCC" 
python optimalT.py

echo "Do surface fitting with the optimal threshold(by MCC) vs dilution rate and coverage"
matlab -nodesktop -nosplash -nodisplay -r "run ./surfacefitting;quit;"

echo Feed the fitted optimal threshold to bayesian test for variant calling
python vcf_Tfit.py
