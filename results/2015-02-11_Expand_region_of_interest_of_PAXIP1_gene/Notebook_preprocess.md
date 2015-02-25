2015-02-11 Expand region of interest of PAXIP1 gene
==============================

Purpose
------------
Find the expanded region for PAXIP1 gene including first three exons and the last exon.  
Rerun rvd on this expanded region to call germline and somatic mutations.

Conclusions
-----------------
Yuting's region: Chromosome 7: 154,738,059 - 154,782,774 (~45kb)  

Expanded region: Chromosome 7:  154,735,400 - 154,794,682 (~59kb)  including three first exons and the last exon.
[http://www.ncbi.nlm.nih.gov/gene?term=NM_007349](http://www.ncbi.nlm.nih.gov/gene?term=NM_007349)


> Why it doesn't include Yuting's region?   
> True region: Chromosome 7: 154,943,687-155,003,084 (~59.40kb)  
> [http://useast.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000157212;r=7:154943687-155003084;redirect=no](http://useast.ensembl.org/Homo_sapiens/Location/View?db=core;g=ENSG00000157212;r=7:154943687-155003084;redirect=no)  


Background
----------------

Materials and Equipment
------------------------------
Average depth for control (HCC1187BL\_S1.dc) is 50x.  
Average depth for case (HCC1187C\_S1) is 70x.

Experimental Protocol
---------------------------
	[fzhang@redwood Expand_region_for_PAXIP1_gene]$ samtools mpileup -C 50 -d 100000  -r chr7:154735400-154794682 -f ../../../../../freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ../../../../../freeze/HCC1187/HCC1187BL_S1.bam > pileup/HCC1187BL_S1.pileup
	[fzhang@redwood Expand_region_for_PAXIP1_gene]$ samtools mpileup -C 50 -d 100000  -r chr7:154735400-154794682 -f ../../../../../freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa ../../../../../freeze/HCC1187/HCC1187C_S1.bam > pileup/HCC1187C_S1.pileup

	[fzhang@redwood bin]$ ./pileup2dc ../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/pileup/HCC1187BL_S1.pileup > ../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/dc/HCC1187BL_S1.dc	
	[fzhang@redwood bin]$ ./pileup2dc ../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/pileup/HCC1187C_S1.pileup > ../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/dc/HCC1187C_S1.dc

	[fzhang@redwood rvd27]/home/pjflaherty/flahertylab/fzhang/Research/rvd/src/python/rvd27
	[fzhang@redwood rvd27]$ python rvd27.py gibbs -o ../../../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/hdf5/HCC1187BL_S1 ../../../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/dc/HCC1187BL_S1.dc -p 10 -s 19860522
	[fzhang@redwood rvd27]$ python rvd27.py gibbs -o ../../../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/hdf5/HCC1187C_S1 ../../../../2015_spring_Fan/RVD2_Review/Expand_region_for_PAXIP1_gene/dc/HCC1187C_S1.dc -p 10 -s 19860522



	
Results
-----------

Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Fan Zhang   Date: 


Witnessed by: ________________________
