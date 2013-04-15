# Program locations here
SAMTOOLS	= /Users/pjflaherty/bin/samtools
PILEUP2DC	= ../../bin/pileup2dc
RVD2		= ../../src/python/rvd27/rvd27.py

# Data folders here
BAMDIR		= ../../data/ms/2012-09-12_Aligned_Data
PILEUPDIR	= ./pileup
DEPTHDIR	= ./depth_chart

# Input Files here
CNTLBAMFILES	= lane1_NoIndex_L001_R1.fastq.gz.prefix.bam \
				  lane2_NoIndex_L002_R1.fastq.gz.prefix.bam
CASEBAMFILES	= lane1_NoIndex_L001_R1.fastq.gz.prefix.bam \
				  lane2_NoIndex_L002_R1.fastq.gz.prefix.bam
				  
FASTQFILE	= ../../data/hg19/hg19.fa.masked
REGION		= chr18:67529192-67625232

# Model output file here
CASEMODELFILE	= case.hdf5

# Nothing below this line needs to change
vpath %.bam $(BAMDIR)
vpath %.pileup $(PILEUPDIR)
vpath %.dc $(DEPTHDIR)

PILEUPFILES	= $(CASEBAMFILES:.fastq.gz.prefix.bam=.pileup)
DEPTHFILES	= $(PILEUPFILES:.pileup=.dc)

run : $(DEPTHFILES)

# Generare the pileup files from BAM files using SAMTOOLS
.INTERMEDIATE : $(PILEUPFILES) $(DEPTHFILES)
$(PILEUPFILES): %.pileup: %.fastq.gz.prefix.bam
	$(SAMTOOLS) mpileup -d 1000000 -r $(REGION) \
		-f $(FASTQFILE) \
		$< \
		> $(addprefix $(PILEUPDIR)/, $@)

# Generate the depth chart files from PILEUP files using pileup2dc
$(DEPTHFILES): %.dc: %.pileup
	$(PILEUP2DC) $(addprefix $(PILEUPDIR)/, $<) > $(addprefix $(DEPTHDIR)/, $@)

# Estimate the model using rvd2.7
.NOTPARALLEL : case
case : $(DEPTHFILES)
	$(RVD2) -vv gibbs -o $(CASEMODELFILE) $(addprefix $(DEPTHDIR)/, $(DEPTHFILES))

# CLeanup the intermediate files
.PHONY : clean
clean :
	rm pileup/*.pileup 
	rm depth_chart/*.dc