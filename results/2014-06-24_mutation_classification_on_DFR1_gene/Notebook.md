2014-06-29 Classification of mutation type in Elisa's Variomics data
==============================

Purpose
------------
To classify the type of mutations in Elisa's Variomics data (Yeast DFR1 gene chrXV:780,906-781,541, length 636). 

Conclusions
-----------------

Background
----------------
For more information about the data please refer to the lab meeting slides from Elisa Lai H.Wong under `\freeze\variomics`.

The reference [DFR1/YOR236W on chromosome XV from coordinates 780906 to 781541](http://www.yeastgenome.org/cgi-bin/getSeq?query=YOR236W&seqtype=ORF%20Genomic%20DNA&format=fasta) is:

	>YOR236W  Chr 15  
	ATGGCTGGAGGAAAGATTCCTATTGTAGGAATTGTGGCATGTTTACAGCCGGAGATGGGG
	ATAGGATTTCGTGGAGGTCTACCATGGAGGTTGCCCAGTGAAATGAAGTATTTCAGACAG
	GTCACTTCATTGACGAAAGATCCAAACAAAAAAAATGCTTTGATAATGGGAAGGAAGACA
	TGGGAATCCATACCGCCCAAGTTTCGCCCACTGCCCAATAGAATGAATGTCATTATATCG
	AGAAGCTTCAAGGACGATTTTGTCCACGATAAAGAGAGATCAATAGTCCAAAGTAATTCA
	TTGGCAAACGCAATAATGAACCTAGAAAGCAATTTTAAGGAGCATCTGGAAAGAATCTAC
	GTGATTGGGGGTGGCGAAGTTTATAGTCAAATCTTCTCCATTACAGATCATTGGCTCATC
	ACGAAAATAAATCCATTAGATAAAAACGCAACTCCTGCAATGGACACTTTCCTTGATGCG
	AAGAAATTGGAAGAAGTATTTAGCGAGCAAGATCCGGCCCAGCTGAAAGAATTTCTTCCC
	CCTAAAGTAGAGTTGCCCGAAACAGACTGTGATCAACGCTACTCGCTGGAAGAAAAAGGT
	TATTGCTTCGAATTCACTCTATACAATCGTAAATGA

We intend to classify deleterious mutations depending on whether they alter the function of essential proteins. Point mutations is small-scale mutations, including:

- **Silent mutations**, which code for the same (or a sufficiently similar) amino acid.
- **Missense mutations**, which code for a different amino acid.
- **Nonsense mutations**, which code for a stop and can truncate the protein.

In the experiment we use [yeast genetic code](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG12) to translate DNA sequence into proteins. 

*[SGD](http://www.yeastgenome.org/) is a good yeast database.*

### Index protocol in different file formats
**One-based index**.Start and end positions are identified using a one-based index. The end position is included. For example, setting start-end to 1-2 describes two bases, the first and second in the sequence.

- BAM files

**Zero-based index**. Start and end positions are identified using a zero-based index. The end position is excluded. For example, setting start-end to 1-2 describes exactly one base, the second base in the sequence.

- BED files
- SNP file

In the experiment, there is a one base forward problem if I use samtools to extract the DFR1 gene (coordinate: chrXV:780906-781541.) from the fasta file  `\freeze\variomics\sacCer2\chrXV.fa` using command
	
	samtools faidx \freeze\variomics\sacCer2\chrXV.fa chrXV:780906-781541 > DFR1.fa

>chrXV:780,906-781,541
CATGGCTGGAGGAAAGATTCCTATTGTAGGAATTGTGGCATGTTTACAGCCGGAGATGGG
GATAGGATTTCGTGGAGGTCTACCATGGAGGTTGCCCAGTGAAATGAAGTATTTCAGACA
GGTCACTTCATTGACGAAAGATCCAAACAAAAAAAATGCTTTGATAATGGGAAGGAAGAC
ATGGGAATCCATACCGCCCAAGTTTCGCCCACTGCCCAATAGAATGAATGTCATTATATC
GAGAAGCTTCAAGGACGATTTTGTCCACGATAAAGAGAGATCAATAGTCCAAAGTAATTC
ATTGGCAAACGCAATAATGAACCTAGAAAGCAATTTTAAGGAGCATCTGGAAAGAATCTA
CGTGATTGGGGGTGGCGAAGTTTATAGTCAAATCTTCTCCATTACAGATCATTGGCTCAT
CACGAAAATAAATCCATTAGATAAAAACGCAACTCCTGCAATGGACACTTTCCTTGATGC
GAAGAAATTGGAAGAAGTATTTAGCGAGCAAGATCCGGCCCAGCTGAAAGAATTTCTTCC
CCCTAAAGTAGAGTTGCCCGAAACAGACTGTGATCAACGCTACTCGCTGGAAGAAAAAGG
TTATTGCTTCGAATTCACTCTATACAATCGTAAATG

This sequence is one base forward comparing with the sequence I obtain from http://www.yeastgenome.org/cgi-bin/getSeq?query=YOR236W&seqtype=ORF%20Genomic%20DNA&format=fasta


Materials and Equipment
------------------------------
- Fasta file  `\freeze\variomics\sacCer2\chrXV.fa`
- VCF files generated in experiment `2014-06-18_apply_rvd27_to_elisa_data`
- mutation classification function `MutClassify.py`
	- The MutClassify works by replace the bases with alternative bases for the whole gene sequence at first, and then access the mutated gene in a whole.
	- This strategy matters when there is positions adjacent mutated together. Is the correct to classify the adjacent point mutations independently, or to classify them together?

Experimental Protocol
---------------------------
1. Generate the fasta file for DFR1 gene using command
	`samtools faidx \freeze\variomics\sacCer2\chrXV.fa chrXV:780907-781542 > DFR1.fa`

	The region of interest is set at chrXV:780907-781542 to avoid the one base forward problem.
2. Please run the program `MutClassify.py` and pass the vcf files generate in experiment `2014-06-18_apply_rvd27_to_elisa_data`


Results
-----------
![](T0S2T1T3POSbar.png)




Future work
------------------------


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 06/24/14


Witnessed by: ________________________
