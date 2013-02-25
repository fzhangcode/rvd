#!/bin/sh
# runall.sh

# Generates a heatmap and fraction genome altered report for TNBC cohort.

# RUNMATLAB: Run the passed matlab function redirecting matlab output. Displays "Failed." or "Success."
runMATLAB () {
	matlab -nodisplay -r "try, eval $1, catch err, exit(1), end, exit(0)" > /tmp/matlab-log-$$.tmp
	return $?
}

printf "Compiling CNV data structure..."
if [ ! -f "CNVData.mat" ]; then runMATLAB "save_cnvdata"; else echo "Skipped."; fi
[ $? -eq 0 ] && echo "Success." || echo " Failed."

printf "Calling CNVs..."
if [ ! -f "CNVCall.mat" ]; then runMATLAB "save_cnvcall"; else echo "Skipped."; fi
[ $? -eq 0 ] && echo "Success." || echo " Failed."

printf "Generating Fraction Genome Altered Table..."
if [ ! -f "fgatable.txt" ]; then runMATLAB "save_fga"; else echo "Skipped."; fi
[ $? -eq 0 ] && echo "Success." || echo " Failed."

printf "Plotting CNV Heatmap..."
if [ ! -f "fgatable.txt" ]; then runMATLAB "save_fga"; else echo "Skipped."; fi
[ $? -eq 0 ] && echo "Success." || echo " Failed."
echo "Display must be active to show and save heatmap."