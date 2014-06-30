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
    ######## Germline mutation: T0 diploid as test sample########
	#CHROM	POS	refBase	altBase	refAminoAcid	altAminoAcid	MutationType
	chrXV	780948	T	.	C	X	Missense Mutation
	chrXV	780959	A	.	E	X	Missense Mutation
	chrXV	780978	T	.	R	X	Missense Mutation
	chrXV	780984	T	.	G	X	Missense Mutation
	chrXV	780997	T	.	L	X	Missense Mutation
	chrXV	781226	A	.	N	X	Missense Mutation
	chrXV	781278	T	.	G	X	Missense Mutation
	chrXV	781403	T	.	V	X	Missense Mutation
	chrXV	781454	T	.	V	X	Missense Mutation
	chrXV	781522	A	.	T	X	Missense Mutation
				
	######## Somatic mutation: T0 diploid as control, T0 haploid as case ########
	#CHROM	POS	refBase	altBase	refAminoAcid	altAminoAcid	MutationType
	chrXV	780921	G	T	K	N	Missense Mutation
	chrXV	781054	A	G	K	E	Missense Mutation
	chrXV	781120	C	A	P	T	Missense Mutation
	chrXV	781246	G	T	E	_	Nonsense Mutation
	chrXV	781348	A	G	K	E	Missense Mutation
	chrXV	781359	T	A	T	T	Silent Mutation
	chrXV	781408	A	G	S	G	Missense Mutation
	chrXV	781450	A	G	K	E	Missense Mutation
	chrXV	781455	A	G	V	V	Silent Mutation
	chrXV	781462	C	G	P	A	Missense Mutation
	chrXV	781518	A	T	E	D	Missense Mutation

	######## Somatic mutation: T0 diploid as control, T1 diploid as case ########
	#CHROM	POS	refBase	altBase	refAminoAcid	altAminoAcid	MutationType
	chrXV	780926	C	T	P	L	Missense Mutation
	chrXV	780950	T	C	L	S	Missense Mutation
	chrXV	780965	G	A	G	E	Missense Mutation
	chrXV	780974	T	C	F	S	Missense Mutation
	chrXV	780991	T	C	W	R	Missense Mutation
	chrXV	781009	A	G	M	A	Missense Mutation
	chrXV	781010	T	C	M	A	Missense Mutation
	chrXV	781012	A	G	K	E	Missense Mutation
	chrXV	781019	T	A	F	Y	Missense Mutation
	chrXV	781021	A	G	R	G	Missense Mutation
	chrXV	781086	A	G	T	T	Silent Mutation
	chrXV	781097	T	C	I	T	Missense Mutation
	chrXV	781108	T	C	F	L	Missense Mutation
	chrXV	781118	T	C	S	P	Missense Mutation
	chrXV	781271	T	G	I	S	Missense Mutation
	chrXV	781292	G	A	S	N	Missense Mutation
	chrXV	781333	A	G	I	V	Missense Mutation
	chrXV	781369	G	A	D	S	Missense Mutation
	chrXV	781370	A	G	D	S	Missense Mutation
	chrXV	781411	G	A	E	K	Missense Mutation
	chrXV	781433	A	C	K	T	Missense Mutation
	chrXV	781457	A	G	E	G	Missense Mutation
	chrXV	781507	T	C	Y	H	Missense Mutation

	######## Somatic mutation: T0 diploid as control, T2 diploid as case ########
	#CHROM	POS	refBase	altBase	refAminoAcid	altAminoAcid	MutationType
	chrXV	780910	G	A	A	T	Missense Mutation
	chrXV	780921	G	A	K	K	Silent Mutation
	chrXV	780931	G	T	V	L	Missense Mutation
	chrXV	780935	G	A	G	E	Missense Mutation
	chrXV	780941	T	C	V	A	Missense Mutation
	chrXV	780944	C	T	A	V	Missense Mutation
	chrXV	780948	T	C	C	C	Silent Mutation
	chrXV	780958	G	A	E	K	Missense Mutation
	chrXV	780972	A	T	G	G	Silent Mutation
	chrXV	780974	T	C	F	S	Missense Mutation
	chrXV	780978	T	G	R	R	Silent Mutation
	chrXV	780984	T	G	G	G	Silent Mutation
	chrXV	780991	T	C	W	R	Missense Mutation
	chrXV	780994	A	C	R	R	Silent Mutation
	chrXV	780997	T	G	L	V	Missense Mutation
	chrXV	781009	A	G	M	A	Missense Mutation
	chrXV	781010	T	C	M	A	Missense Mutation
	chrXV	781012	A	G	K	E	Missense Mutation
	chrXV	781019	T	C	F	S	Missense Mutation
	chrXV	781097	T	C	I	T	Missense Mutation
	chrXV	781108	T	C	F	L	Missense Mutation
	chrXV	781118	T	C	S	P	Missense Mutation
	chrXV	781136	T	A	V	D	Missense Mutation
	chrXV	781200	T	C	S	S	Silent Mutation
	chrXV	781204	T	G	S	A	Missense Mutation
	chrXV	781226	A	T	N	I	Missense Mutation
	chrXV	781268	T	A	V	E	Missense Mutation
	chrXV	781278	T	C	G	G	Silent Mutation
	chrXV	781333	A	G	I	V	Missense Mutation
	chrXV	781335	A	G	I	V	Missense Mutation
	chrXV	781369	G	A	D	N	Missense Mutation
	chrXV	781403	T	C	V	A	Missense Mutation
	chrXV	781406	T	G	F	C	Missense Mutation
	chrXV	781454	T	A	V	E	Missense Mutation
	chrXV	781468	A	T	T	S	Missense Mutation
	chrXV	781482	A	T	Q	H	Missense Mutation
	chrXV	781506	T	C	G	G	Silent Mutation
	chrXV	781522	A	T	T	S	Missense Mutation
	chrXV	781529	A	T	Y	F	Missense Mutation

Future work
------------------------


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 06/24/14


Witnessed by: ________________________
