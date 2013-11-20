

PILEUP2DC=../../bin/pileup2dc

BAMDIR := ../2013-10-11_new_dataset_searching/HCC1187
BAMFILE := $(wildcard $(addprefix $(BAMDIR)/, *.bam)) # all pair1 & pair2

PILEUPDIR = ./pileup
PILEUPFILE := $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(notdir $(BAMFILE) ) ) ) )
REGION = chr7:154738059-154782774
FASTAFILE = ../../data/freeze/hg19/hg19.fa.masked

DEPTHDIR = ./depth_chart
DEPTHFILE := $(addprefix $(DEPTHDIR)/, $(addsuffix .dc, $(basename $(notdir $(PILEUPFILE) ) ) ) )

.INTERMEDIATE: %.sorted.bam %.pileup

all: $(DEPTHFILE)
	echo $(BAMFILE)

clean: 
	rm $(PILEUPFILE)

#TODO: form correct prerequisite

# Create the pileup file
$(PILEUPFILE): $(BAMFILE)
	mkdir -p $(PILEUPDIR)
	samtools mpileup \
		-d 1000000 \
		-r $(REGION) \
	        -f $(FASTAFILE) \
		$(addprefix $(BAMDIR)/, $(patsubst %.pileup, %.bam, $(@F)) ) \
		> $@

# Create the depth chart file
$(DEPTHFILE): $(PILEUPFILE)
	mkdir -p $(DEPTHDIR)
	$(PILEUP2DC) $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(@F)))) > $@
