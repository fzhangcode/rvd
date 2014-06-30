2014-06-30 Output the variants information into correctly formatted vcf files. 
==============================

Purpose
------------
Output the variants called by rvd27 into correctly formatted vcf files. 

Conclusion
-----------
IGV is able to load vcf files we generated from rvd27.

Background
-----------
Please find the manual for **VCF (Variant Call Format) version 4.0** at [http://www.1000genomes.org/node/101](http://www.1000genomes.org/node/101 "Manual")

Experimental Protocol
-----------
Load variant calls (.vcf files) into IGV
>1. Index vcf file
	1. Choose Tools ? Run igvtools....
	2. Choose "index" from the commands drop-down menu
	3. Select your *.vcf file (Ex: filename.vcf) for "Input File"
	4. Click the "run" button.
	5. Click the "close" button to close the panel. A "filename.vcf.idx" file has generated at the vcf file directory
>2. Load vcf file
	1. choose File-> Load from File. 
	2. Navigate to your vcf file directory and load your indexed vcf file. 

*Tip: You can also index BAM and FASTA files the same way inside of IGV if you haven't already created indexes for them. But, it's usually easier and quicker to do this on the command line.*

Results
-----------
IGV is able to load vcf files we generated from rvd27.

Future work
------------------------

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 06/27/14


Witnessed by: ________________________
