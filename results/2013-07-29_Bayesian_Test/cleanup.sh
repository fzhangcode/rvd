#!/bin/sh
# cleanup.sh

# Script to clean up temporary data files so results can be regenerated.

printf "Cleaning up temporary files and directories..."
rm -rf *.hdf5 *.pdf
printf "Done.\n"

exit 0
