2013-10-01 SNP calling using strelka
==============================

Purpose
------------
To compare strelka with rvd2 on variant calling

Conclusions
-----------------
Working experiment.

Background
----------------
Strelka is an analysis package designed to detect somatic SNVs and small indels from the aligned sequencing reads of matched tumor-normal samples. Webpage
[Strelka somatic variant caller - Google Sites](https://sites.google.com/site/strelkasomaticvariantcaller/) contains all the information about Strelka.


Materials and Equipment
------------------------------



Experimental Protocol
---------------------------
To avoid the trouble caused by a redundant line in bam files header, please run `header_convert.sh` in `2013-10-02_problematic_header_removal` to generate the reheadered bam files first.

Run `SNP_strelka.sh` in current directory to get the vcf files, which will be available under directory 'work'.


Results
-----------



Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: _____Yuting He____________     Date: _______10/02/13______________


Witnessed by: ________________________
