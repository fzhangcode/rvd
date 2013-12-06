
#PICARD=/usr/local/pjf/picard-tools-1.96/DownsampleSam.jar
#PICARD=java -Xmx2g -jar ../2013-08-06_Downsample_Read_Depth/picard-tools-1.96/DownsampleSam.jar

PILEUP2DC=../../bin/pileup2dc

BAMDIR := ../../data/HCC1187
BAMFILE := $(wildcard $(addprefix $(BAMDIR)/, *.bam)) # all pair1 & pair2

PILEUPDIR = ./pileup
PILEUPFILE := $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(notdir $(BAMFILE) ) ) ) )
REGION = chr7:154738059-154782774
FASTAFILE = ../../data/freeze/hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

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
