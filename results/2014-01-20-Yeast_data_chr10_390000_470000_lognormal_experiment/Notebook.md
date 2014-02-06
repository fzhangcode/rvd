2014-01-14 Test RVD3 with log-normal prior on the yeast data (Chr10)
==============================

Purpose
------------
Call variants between Chr10 390000:470000

Call the gene CYR1. CYR1 is found in Chr10 on position 428286 (Ref:A->Alt:C) by kvitek2013 (supplement Table_S2).

Conclusions
-----------------
The log-normal prior on RVD3 works well because in calls mutations without false positive.

Background
-----------------


Materials and Equipment
------------------------------
Run original data in freeze/kvitek2013. The generated hdf5 and vcf files are saved in the current directory.

Experimental Protocol
---------------------------
Run runnall.py to call rvd29.py

Run test_rvd29.py 

Results
-----------
From the VCF file generated, we call 27 mutations without false positive. And all these variants are covered in the results from Jeffreys prior except one (Position=453972)

But we failed to find the CYR1 because it doesn't show in the dc file neither.

Archived Samples
-------------------------
All the original data are stored in folder: freeze/kvitek2013

Archived Computer Data
------------------------------


Prepared by: ____  Fan Zhang________     Date: ______2014/01/24  __


Witnessed by: ________________________
