2013-10-02 remove problematic line in header
==============================

Purpose
------------
To remove a line in header which causes problem to Virmid and strelka

Conclusions
-----------------
Working code available in the same directory as `header_convert.sh`

Background
----------------
The header of the bam files contain:

	@HD     VN:1.4  SO:coordinate
	@SQ     SN:1_1_400      LN:400
	@SQ     SN:2_1_10       LN:10
	@PG     ID:bwa  PN:bwa  VN:0.5.9-r16

The second 'ref' `2_1_10`  was said to be a dummy required by the read depth chart script (from John Bell) and no longer needed. Therefore, it was directly removed from the bam files.

Materials and Equipment
------------------------------
`header_convert.sh`


Experimental Protocol
---------------------------

Results
-----------


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____Yuting He____________     Date: _______10/02/13______________


Witnessed by: ________________________
