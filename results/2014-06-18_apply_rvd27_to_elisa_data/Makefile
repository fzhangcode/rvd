PILEUP2DC=../../bin/pileup2dc
PILEUP2DCMAIN=../../src/pileup2dc/main.c

### specify your bamfile directory here #add directory and file name
BAMFILE = ../../data/variomics/bamfolder/T0diploid_S2.bam
FASTAFILE = ../../data/variomics/sacCer2/chrXV.fa

## options for wildcard
# # BAMDIR = ../../data/GSY11351
# # BAMFILE := $(wildcard $(addprefix $(BAMDIR)/, *.bam)) 
# # FASTAFILE = $(addprefix $(BAMDIR)/, plasmid.fa)

SORTEDBAM := $(patsubst %.bam, %.sorted.bam, $(BAMFILE) )
INDEXBAM = $(addsuffix .bai, $(SORTEDBAM))

PILEUPDIR = ./pileup
PILEUPFILE := $(addprefix $(PILEUPDIR)/, $(addsuffix .pileup, $(basename $(notdir $(BAMFILE) ) ) ) )
REGION = chrXV:780,906-781,541 


DEPTHDIR = ./depth_chart/
DEPTHFILE := $(addprefix $(DEPTHDIR)/, $(addsuffix .dc, $(basename $(notdir $(PILEUPFILE) ) ) ) )

all: $(DEPTHFILE)

clean: 
	rm $(SORTEDBAM)
	rm $(PILEUPFILE)
	rm $(DEPTHFILE)


## Need to make a loop if there is wild card?
# sort the BAM file
$(SORTEDBAM): $(BAMFILE)
	samtools sort $^ $(basename $@)

# index sorted bam file
$(INDEXBAM): $(SORTEDBAM)
	samtools index $^ > $@

# Create the pileup file
$(PILEUPFILE): $(SORTEDBAM) $(INDEXBAM)
	mkdir -p $(PILEUPDIR)
	samtools mpileup -C 50 -d 1000000 -f $(FASTAFILE) -r $(REGION) $(SORTEDBAM) > $@

$(PILEUP2DC): $(PILEUP2DCMAIN)
	gcc -o $@ $^

# Create the depth chart file
$(DEPTHFILE): $(PILEUP2DC) $(PILEUPFILE)
	mkdir -p $(DEPTHDIR)
	$^ > $@

