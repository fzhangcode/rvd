# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
from __future__ import print_function
from __future__ import division
import numpy as np
import pdb

import matplotlib.pyplot as plt
import pandas as pd
import vcf
import re

import random
import textwrap
import os

def main():
    ## read in the fasta files
    filename = './../2014-06-18_apply_rvd27_to_elisa_data/DFR11.fa'
    with open(filename,'r') as fin:
        RefSeq = fin.read()
    RefSeq = clean_sequence(RefSeq)
    RefPepSeq= translate(RefSeq)

    ROI = range(780907,781543)

    ######## Germline mutation: T0 diploid as test sample########
    vcffile = './../2014-06-18_apply_rvd27_to_elisa_data/T0diploid_S2_0_01.vcf'
    [AltPepSeq, pos, refb, altb] = import_altSeq(vcffile, RefSeq, ROI, srows = 9)
    mutationclf(RefPepSeq, AltPepSeq, pos, ROI, refb, altb, outputfile = os.path.basename(vcffile).replace('vcf','txt'))

    # ######## Somatic mutation: T0 diploid as control, T0 haploid as case ########
    # vcffile = './../2014-06-18_apply_rvd27_to_elisa_data/T0haploid_S1.vcf'
    # [AltPepSeq, pos, refb, altb] = import_altSeq(vcffile, RefSeq, ROI)
    # mutationclf(RefPepSeq, AltPepSeq, pos, ROI, refb, altb, outputfile = os.path.basename(vcffile).replace('vcf','txt'))

    # ######## Somatic mutation: T0 diploid as control, T1 diploid as case ########
    # vcffile = './../2014-06-18_apply_rvd27_to_elisa_data/T1diploid_S3.vcf'
    # [AltPepSeq, pos, refb, altb] = import_altSeq(vcffile, RefSeq, ROI)
    # mutationclf(RefPepSeq, AltPepSeq, pos, ROI, refb, altb, outputfile = os.path.basename(vcffile).replace('vcf','txt'))

    # ######## Somatic mutation: T0 diploid as control, T2 diploid as case ########
    # vcffile = './../2014-06-18_apply_rvd27_to_elisa_data/T2diploid_S4.vcf'
    # [AltPepSeq, pos, refb, altb] = import_altSeq(vcffile, RefSeq, ROI)
    # mutationclf(RefPepSeq, AltPepSeq, pos, ROI, refb, altb, outputfile = os.path.basename(vcffile).replace('vcf','txt'))    

    ## what about bases next to each other?
def import_altSeq(vcffile, RefSeq, ROI, srows = 10):
## read in the vcf files
    
    t = pd.read_table(vcffile, header = 1, skiprows = srows)
    pos = t.POS
    pos = [int(s) for s in pos]

    altb = t.ALT
    refb = t.REF

    ## create the mutated sequence
    AltSeq = list(RefSeq)
    # pdb.set_trace()

    for i in pos:
        AltSeq[ROI.index(i)] = altb[pos.index(i)]

    ## convert the AltSeq to a str
    AltSeq = ''.join(s for s in AltSeq)
    AltPepSeq = translate(AltSeq)
    return AltPepSeq, pos, refb, altb


def mutationclf(sequence1, sequence2, pos, ROI, refb, altb, outputfile = 'Mutationtypes.txt'):
    f = open(outputfile, 'w')
    print('##The original peptide sequence is:', file = f)
    # print('%s' %sequence1, file = f)
    f.write(textwrap.fill(sequence1, width = 50))
    print('\n##The mutated peptide sequence is:', file = f)
    # print('%s' %sequence2, file = f)
    f.write(textwrap.fill(sequence2, width = 50))

    sequence1 = list(sequence1)
    sequence2 = list(sequence2)
   
    print('\n#CHROM\tPOS\trefBase\taltBase\trefAminoAcid\taltAminoAcid\tMutationType', file = f)

    for i in pos:
        pep_idx = int(ROI.index(i)/3)
        pos_idx = pos.index(i)
        # pdb.set_trace()
        if sequence1[pep_idx] == sequence2[pep_idx]:
            print('chrXV\t%d\t%s\t%s\t%s\t%s\tSilent Mutation' %(i, refb[pos_idx], altb[pos_idx], sequence1[pep_idx], sequence2[pep_idx]), file = f)
        elif sequence1[pep_idx] != sequence2[pep_idx] and sequence2[pep_idx]!= '_':
            print('chrXV\t%d\t%s\t%s\t%s\t%s\tMissense Mutation' %(i, refb[pos_idx], altb[pos_idx], sequence1[pep_idx], sequence2[pep_idx]), file = f)
        elif sequence2[pep_idx]== '_': 
            print('chrXV\t%d\t%s\t%s\t%s\t%s\tNonsense Mutation' %(i, refb[pos_idx], altb[pos_idx], sequence1[pep_idx], sequence2[pep_idx]), file = f)

    f. close()


def random_sequence( length ):
    """Return a random string of AGCT of size length."""
    return ''.join([ random.choice( 'AGCT' ) for i in range(int( length ))])

def clean_sequence( sequence ):
    """Given a sequence string, return a crap-free, standardized DNA version."""
    s = sequence.replace( '\r', '' ).split( '\n' )  # separate each line
    if s[0][0] == '>': s = s[ 1 :]                  # remove defline
    s = ''.join( s )                                # make one long string
    s = s.replace( ' ', '' ).replace( '\t', '' )    # remove spaces
    return s.upper().replace( 'U', 'T' )

def report_bad_chars( sequence ):
    """Given a string 'sequence', return a dictionary of any non-AGCT characters."""
    bad_chars = {}
    for l in sequence:
        if l not in 'AGCT':
            if l in bad_chars: bad_chars[ l ] += 1
            else: bad_chars[ l ] = 1
    if bad_chars != {}: print( bad_chars )


def translate( sequence ):
    """Return the translated protein from 'sequence' assuming +1 reading frame"""
    ## Yeast genetic code
    # Difference from the standard Code: CUG-> Serine; standard CUG->Leucine
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'S', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    return ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])

if __name__ == "__main__":
    main()