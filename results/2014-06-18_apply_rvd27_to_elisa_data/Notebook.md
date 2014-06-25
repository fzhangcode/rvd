2014-06-18 Apply RVD27 to Elisa's Variomics data
==============================

Purpose
------------
To detect variants in Elisa's Variomics data (Yeast DFR1 gene chrXV:780,906-781,541, length 635). 

Conclusions
-----------------

Background
----------------
Please refer to the lab meeting slides from Elisa Lai H.Wong under `\freeze\variomics`.

Materials and Equipment
------------------------------
Reference fasta file was download from [here ](http://hgdownload-test.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/ ). After decompression, file `chrXV.fa` was used as reference file.

Bam files under `\freeze\variomics\bamfolder` were processed. The files were download from [here](http://chemogenomics.pharmacy.ubc.ca/MAU/NGS_LINKS/A9JEK/A9JEK/).

	username: twinkle
	password: littlestar
More data details please refer to the Notebook under `\freeze\variomics`


Experimental Protocol
---------------------------
1. Run `2014-06-18_apply_rvd27_to_elisa_data\Makefile` to generate depth chart files from bam files.

	parameter setting:

		fasta file: \freeze\variomics\sacCer2\chrXV.fa
		Region of interest: chrXV:780,906-781,541 
		In samtools mpileup: -C 50 -- Alignment using BWA; -d 1,000,000 maximum reads per input BAM
	You need to specify the bam file directory and file name in `BAMFILE = ../../data/variomics/bamfolder/XXX.bam` to run each bamfile.
2. Run bash file `2014-06-18_apply_rvd27_to_elisa_data\hdf5gen.sh` to generate hdf5 files using the `gibbs` function in `rvd27` program. 

	parameter setting:

		gibbs sampling size: 4000; metropolis-hastings sampling size 5
3. Run the variant detection functions in `rvd27` program
	The bash files `variant_somatic_test.sh` and `variant_germline_test.sh` are bash files which can be run given correct control files and case file.

	Germline test on `T0diploid_S2.bam` sample. 

		There were three detection MAF threshold set in the experiment: 0.001, 0.005, 0.01
		The proper threshold should be:
		The credible level alpha was set at 0.05.

	There were four sets of somatic test to test if there are positions that have significantly different MAF level, including

		control sample:				case sample					output vcf files
		T0diploid_S2.bam			T0haploid_S1.bam			T0haploid_S1.vcf 
		T0diploid_S2.bam			T1diploid_S3.bam			T1diploid_S3.vcf	
		T0diploid_S2.bam			T2diploid_S4.bam			T2diploid_S4.vcf	
		T1diploid_S3.bam			T2diploid_S4.bam			T1diploid_S3_T2diploid_S4.vcf

	The detection threshold tau was set at 0.0, and the credible level alpha was set at 0.05. 

Results
-----------
![](T0S2T1T3POSbar.png)

**Figure: Bar plots shows the MAF level of positions which has significantly different MAF levels in T2 generation and T1 generation from T0 (diploid) generation. There were 11 positions reported, positions in common in file `T1diploid_S3.vcf` and file `T2diploid_S4.vcf`.**

A general trend in the bar plot shows that the MAF levels in these 11 positions first increase (from generation T0 to T1), but then decrease to the similar MAF level as generation T0 (from generation T1 to T2). 

![](T1_T0S2.png)

**Figure: MAF levels of positions where MAF levels differ significantly in T1 generation from T0 generation (diploid). There are 23 positions reported, as shown in file `T1diploid_S3.vcf`**

It can be seen that the MAF levels increased in T1 generation compared with T0 generation in all positions. 

![](T2_T0S2.png)

**Figure: MAF levels of positions where MAF levels differ significantly in T2 generation from T0 generation (diploid). There are 39 positions reported, as shown in file `T2diploid_S4.vcf`**

It can be seen that the MAF levels did not increase universally across all the reported positions. MAF level in 25 positions decreased and increased in the remaining 14 positions. The decreasing effect is more noticeable than the increasing effect. 


![](T2_T1.png)

**Figure: MAF levels of positions where MAF levels differ significantly in T2 generation from T1 generation. There are 17 positions reported.**

It can be seen that there is a deceasing trend in the MAF levels from T1 generation to T2 generation. 

**Elisa diluted the pool after T1(20) generation, this might explain why the MAF level decreases affterwards.**

![](vennfig.png)

**Figure: Venn diagrams for variants detected using different tests and samples. The different sets in the venn diagrams were labels in the format as `CaseSample-ControlSample` for somatic mutations and `ControlSample` for germline mutations.**



![](T0S2_T0S1.png)

**Figure: MAF levels of positions where MAF levels differ significantly in T0 generation haploid sample from T0 generation diploid sample. There are 11 positions reported.**

It can be seen that the MAF levels in haploid sample is universally higher than the MAF levels in the diploid sample. 

### Bioinformatics Analysis
		control sample:				case sample					output vcf files
		T0diploid_S2.bam			T0haploid_S1.bam			T0haploid_S1.vcf 
		T0diploid_S2.bam			T1diploid_S3.bam			T1diploid_S3.vcf	


Future work
------------------------
DFR1 region with extension, Chromosome 12, from base 779906 to base 782541 (length 2635) with no upstream or downstream probe.

Translate DNA to Protein

[sequence for a region of DFR1/YOR236W](http://www.yeastgenome.org/cgi-bin/getSeq?query=YOR236W&seqtype=ORF%20Genomic%20DNA&format=fasta)


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

deleterious mutations

[yeast genetic code](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c)

start/stop codon
yeast database
http://www.yeastgenome.org/

A **nonsense mutation** is a point mutation in a sequence of DNA that results in a premature stop codon, or a nonsense codon in the transcribed mRNA, and possibly a truncated, and often nonfunctional protein product. (See Stop codon.)

**Missense mutations** or nonsynonymous mutations are types of point mutations where a single nucleotide is changed to cause substitution of a different amino acid. This in turn can render the resulting protein nonfunctional. Such mutations are responsible for diseases such as Epidermolysis bullosa, sickle-cell disease, and SOD1 mediated ALS (Boillée 2006, p. 39).

A **neutral mutation** is a mutation that occurs in an amino acid codon that results in the use of a different, but chemically similar, amino acid. The similarity between the two is enough that little or no change is often rendered in the protein. For example, a change from AAA to AGA will encode arginine, a chemically similar molecule to the intended lysine.

**Silent mutations** are mutations that do not result in a change to the amino acid sequence of a protein, unless the changed amino acid is sufficiently similar to the original. They may occur in a region that does not code for a protein, or they may occur within a codon in a manner that does not alter the final amino acid sequence. The phrase silent mutation is often used interchangeably with the phrase synonymous mutation; however, synonymous mutations are a subcategory of the former, occurring only within exons (and necessarily exactly preserving the amino acid sequence of the protein). Synonymous mutations occur due to the degenerate nature of the genetic code


Archived Samples
-------------------------

Archived Computer Data
------------------------------


Prepared by: Yuting He     Date: 06/18/14


Witnessed by: ________________________
