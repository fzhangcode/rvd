#!/bin/sh

for f in `find ./Newfolder/ -type f`; do
	#statements
	mv -v $f $f.vcf
done
