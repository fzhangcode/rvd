
PICARD=java -Xmx2g -jar /usr/local/pjf/picard-tools-1.96/DownsampleSam.jar
PILEUP2DC=../../bin/pileup2dc

BAMDIR := /research/pjflaherty/flahertylab/freeze/h1n1/Synthetic_BAM_files
BAMFILE := $(wildcard $(addprefix $(BAMDIR)/, 20100916_c?_p?.*.bam)) # all pair1 & pair2

DRATE = 0.000333333
DFRAC = 3000
DOWNSAMPLEDIR = ./bam/$(DFRAC)
DOWNSAMPLEBAM := $(addprefix $(DOWNSAMPLEDIR)/, $(notdir $(BAMFILE) ) )
# := enforces wildcard expansion

DOWNSAMPLESORTEDBAM := $(patsubst %.bam, %.sorted.bam, $(DOWNSAMPLEBAM) )

PILEUPDIR = ./pileup/$(DFRAC)
PILEUPFILE := $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(notdir $(DOWNSAMPLEBAM) ) ) ) )
REGION = chr1:1-400
FASTAFILE =  /research/pjflaherty/flahertylab/freeze/h1n1/Synthetic_BAM_files/plasmid.fa

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
