
PICARD=java -Xmx2g -jar /usr/local/pjf/picard-tools-1.96/DownsampleSam.jar
#PICARD=java -Xmx2g -jar ../2013-08-06_Downsample_Read_Depth/picard-tools-1.96/DownsampleSam.jar
PILEUP2DC=../../bin/pileup2dc

BAMDIR := ../../../../../freeze/baker_yeast/GSY1135/test/right_bam/
BAMFILE := $(wildcard $(addprefix $(BAMDIR)/, *.bam)) # all pair1 & pair2

DRATE = 0.1
DFRAC = 10
DOWNSAMPLEDIR = ./bam/$(DFRAC)
DOWNSAMPLEBAM := $(addprefix $(DOWNSAMPLEDIR)/, $(notdir $(BAMFILE) ) )
# := enforces wildcard expansion

DOWNSAMPLESORTEDBAM := $(patsubst %.bam, %.sorted.bam, $(DOWNSAMPLEBAM) )

PILEUPDIR = ./pileup/$(DFRAC)
PILEUPFILE := $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(notdir $(DOWNSAMPLEBAM) ) ) ) )
REGION = chr1:1-400
FASTAFILE = ../../../../../freeze/baker_yeast/GSY1135/Reference_GSY1135/GSY1135_Chr10.fasta

DEPTHDIR = ./depth_chart/$(DFRAC)
DEPTHFILE := $(addprefix $(DEPTHDIR)/, $(addsuffix .dc, $(basename $(notdir $(PILEUPFILE) ) ) ) )

all: $(DEPTHFILE)

clean: 
	rm $(DOWNSAMPLEBAM)
	rm $(DOWNSAMPLESORTEDBAM)
	rm $(PILEUPFILE)

#TODO: form correct prerequisite
# downsample the BAM file
$(DOWNSAMPLEBAM): $(BAMFILE)
	mkdir -p $(DOWNSAMPLEDIR)
	$(PICARD) \
		INPUT=$(addprefix $(BAMDIR)/, $(@F) ) \
		OUTPUT=$@ \
		RANDOM_SEED=10241978 \
		PROBABILITY=$(DRATE) \
		VALIDATION_STRINGENCY=SILENT

# sort the downsampled BAM file
$(DOWNSAMPLESORTEDBAM): $(DOWNSAMPLEBAM)
	samtools sort $(patsubst %.sorted.bam, %.bam, $@) $(basename $@)

# Create the pileup file
$(PILEUPFILE): $(DOWNSAMPLESORTEDBAM)
	mkdir -p $(PILEUPDIR)
	samtools mpileup \
                -C 50 \
		-d 1000000 \
		-f $(FASTAFILE) \
		$(addprefix $(DOWNSAMPLEDIR)/, $(patsubst %.pileup, %.sorted.bam, $(@F)) ) \
		> $@

# Create the depth chart file
$(DEPTHFILE): $(PILEUPFILE)
	mkdir -p $(DEPTHDIR)
	$(PILEUP2DC) $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(@F)))) > $@
