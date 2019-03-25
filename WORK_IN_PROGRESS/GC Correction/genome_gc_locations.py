#!/usr/bin/python

'''
Takes a bed file as input and outputs the a bed file with GC content 
of regions added.
'''

import sys
import os
import decimal
import argparse
from Bio import SeqIO


def gc_content(seq, seqlen):
	return decimal.Decimal(seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')) / seqlen

parser = argparse.ArgumentParser(description = "Takes a bed file as input and outputs the a bed file with GC content of regions added")
parser.add_argument("input", help = "BED file to add GC content to")
parser.add_argument("reference", help = "FASTA file of reference genome of BED file")
parser.add_argument("output", help = "Output filename")
parser.add_argument("-b", "--basepairs", help = "Size of genomic window (in bp) to calculate GC content (default = 500)", type = int, default = 500)
parser.add_argument("-g", "--gcwindow", help = "GC content bin size (default = 0.05)", default = '0.05')

args = parser.parse_args()

gc_winsize = decimal.Decimal(args.gcwindow)
gc_places = gc_winsize.as_tuple().exponent
genome_winsize = args.basepairs

infile = open(args.input, 'r')
reference = SeqIO.index(args.reference, 'fasta')
outfile = open('%s_TEMPFILE' %(args.output), 'w')

gc_loc_counts = [0] * (int(1/gc_winsize) + 1)

current_chrom = None
current_gc = -1

for line in infile:
	if line.startswith('#'):
		continue
	data = line.rstrip().split()
	if data[0] != current_chrom:
		current_chrom = data[0]
		reference_chrom = reference[current_chrom]
	for base in xrange(int(data[1]), int(data[2])):
		gc = gc_content(reference_chrom.seq[base: base + genome_winsize], genome_winsize)
		gc = int(gc * (1/gc_winsize)) * gc_winsize
		gc_win =  int(gc/gc_winsize)
		gc_loc_counts[gc_win] += 1
		if current_gc != gc:
			if current_gc >= 0:
				outfile.write('%s\t%i\t%i\t%.*f\n' %(current_chrom, start_base, base, -gc_places, current_gc))
			current_gc = gc
			start_base = base
	outfile.write('%s\t%i\t%i\t%.*f\n' %(current_chrom, start_base, base+1, -gc_places, gc))
	current_gc = -1

infile.close()
reference.close()
outfile.close()
outfile = open(args.output, 'w')
outfile.write('#Genome Window Size: %i bp\tGC Content Window Size: %.*f\n' %(genome_winsize, -gc_places, gc_winsize))
for idx, window in enumerate(gc_loc_counts):
	outfile.write('#GC Content of %.*f-%.*f: %i locations\n' %(-gc_places + 2, idx * gc_winsize, -gc_places + 2, min(1., (idx + 1) * gc_winsize - decimal.Decimal((0, (0,1), gc_places - 2))), window))
tempout = open('%s_TEMPFILE' %(args.output), 'r')
for line in tempout:
	outfile.write(line)

outfile.close()
tempout.close()
os.remove('%s_TEMPFILE' %(sys.argv[3]))
