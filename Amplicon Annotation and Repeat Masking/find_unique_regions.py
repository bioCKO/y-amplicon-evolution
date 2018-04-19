#!/usr/bin/python

'''
Takes FASTA file of a genome or chromosome as input and finds the uniquely mappable regions.

This script works by generating all possible reads of a given size from a reference genome 
or chromosome by moving a window of that size one base at a time over each chromosome.
Those reads are then mapped to the entire reference genome using bowtie2.
Finally, any reads which either don't map anywhere (because of the presence of Ns in the reads) 
or that map to more than one place with an chosen alignment score or higher are detected, 
and the program returns a mapability file in the bed format.

NOTE: This script takes a very long time to run.
'''

#TO DO: Add option for alternatives to using Bowtie2

from Bio import SeqIO
import subprocess
import argparse
import datetime
import os

parser = argparse.ArgumentParser(description = "Takes FASTA file of a genome or chromosome as input and returns a mappability file of the uniquely mappable regions")
parser.add_argument("reference_genome", help = "FASTA file of the genome or chromosome")
parser.add_argument("window_size", help = "Size of window to check mapability; set this equal to the read length of your sequencing data") 
parser.add_argument("reference_index", help = "Bowtie2 index prefix of the genome")
parser.add_argument("-d", "--output_directory", default = "./", help = "Directory to output the mapability file")
parser.add_argument("-m", "--min_score", type = int, default = -14, help = "Minimum alignment score to treat as multimapping (default: -14)")
parser.add_argument("-o", "--output", help = "Custom output filename")

args = parser.parse_args()

out_dir = args.output_directory if args.output_directory.endswith('/') else args.output_directory + '/'
if args.output:
	mapname = args.output
else:
	mapname = out_dir + '.'.join(args.reference_genome.split('/')[-1].split('.')[:-1]) + '_window_' + args.window_size + '_mapability.bed' if '.' in args.reference_genome and args.reference_genome.split('.')[-1] in ['fasta', 'fa'] else out_dir + args.reference_genome.split('/')[-1] + '_window_' + args.window_size + '_mapability.bed'
	#The default name of the final mapability file: if name of input file is x.fa or x.fasta, will be x_window_{window_size}_mapability.bed, otherwise will be full input filename followed by _window_{window_size}_mapability.bed
if os.path.isfile(mapname):
	print 'The output file name already exists. Remove that file, or use a different directory.'
	quit()

mapfile = open(mapname, 'w')

ref_genome = SeqIO.parse(args.reference_genome, 'fasta')
for chromosome in ref_genome:
	#Set names of temporary files: the file of all possible reads and the sam file of those reads aligned to the reference genome
	read_out = out_dir + chromosome.id + '_'.join(str(datetime.datetime.now()).split()) + '_all_possible_reads.temp'
	sam_out =  out_dir + chromosome.id + '_'.join(str(datetime.datetime.now()).split()) + '_all_possible_reads_aligned.sam.temp'
	
	#Step 1: Generate all possible reads from a chromosome

	print 'Generating all possible reads from %s...' %(chromosome.id)
	with open(read_out, 'w') as outfile:
		readcount = 1
		for read in xrange(len(chromosome) - int(args.window_size) + 1):
			outfile.write('>%sr%i\n' %(chromosome.id, readcount))
			outfile.write(str(chromosome.seq[read:read + int(args.window_size)]) + '\n')
			readcount += 1

	#Step 2: Align the generated reads back to the whole reference genome

	print 'Aligning %s reads to reference...' %(chromosome.id)
	print subprocess.check_output('bowtie2 -f -k 2 --score-min C,%i -x %s -U %s -S %s' %(args.min_score, args.reference_index, read_out, sam_out), shell = True),
	print subprocess.check_output('rm -f %s' %(read_out), shell = True),
	
	#Step 3: Detect uniquely mappable regions of the genome from the alignment
	
	print 'Finding uniquely mappable regions for %s...' %(chromosome.id)
	with open(sam_out, 'r') as samfile:
		start = 1 #The starting base for the current uniquely mappable region
		for read in samfile:
			#Skip header lines of SAM file
			if read.startswith('@'):
				continue
			samread = read.split()
			window_position = int(samread[0].split('r')[-1])
			#Check if read is unmapped or mapped to more than one location using the binary SAM flag; if it is, end the current uniquely mappable region and write it to the mapability file
			if int(samread[1]) & 260 != 0:
				if window_position != start: #Check if the previous base was also unmappable; if so, the mappable region ended by this unmappable read has length 0, so skip it
					mapfile.write('%s\t%i\t%i\n' %(chromosome.id, start - 1, window_position - 1))
				start = window_position + 1 #Set the start of a new mappable region to the next base
		if window_position != start - 1: #After going over the entire chromosome, check if the final bases were uniquely mappable, and if so write them to the mapability file
			mapfile.write('%s\t%i\t%i\n' %(chromosome.id, start - 1, window_position))
	print subprocess.check_output('rm -f %s' %(sam_out), shell = True),
mapfile.close()