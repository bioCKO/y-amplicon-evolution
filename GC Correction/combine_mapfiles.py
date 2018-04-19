#!/usr/bin/python

'''
Combine GC-content mapability files into one.
'''

import argparse
import decimal
import os

parser = argparse.ArgumentParser(description = '')
parser.add_argument('infiles', nargs='*', help = 'GC-added mapability files to combine')
parser.add_argument('-o', '--output', default = 'genome_mapability.bed', help = 'Output filename (default = genome_mapability.bed)')
parser.add_argument('-w', '--override_window_mismatch', action = 'store_true', help = 'Allow combining of files with different genome window sizes')

args = parser.parse_args()

if os.path.isfile(args.output):
	print 'Output file already exists. Please choose another ouput filename.'
	quit()
	

if len(args.infiles) < 2:
	print 'Enter at least two files to combine'
	quit()

for fileno, infilename in enumerate(args.infiles):
	infile = open(infilename, 'r')
	gc_linecount = 0
	for line in infile:
		if line.startswith('#Genome Window Size'):
			file_genome_winsize = int(line.split()[3])
			file_gc_winsize = decimal.Decimal(line.split()[-1])
			if fileno == 0:
				genome_winsize = file_genome_winsize
				gc_winsize = file_gc_winsize
				gc_places = gc_winsize.as_tuple().exponent
				gc_loc_counts = [0] * (int(1/gc_winsize) + 1)
			else:
				if genome_winsize != file_genome_winsize:
					print 'ERROR: Genome Window Sizes of files do not match!'
					if not args.override_window_mismatch:
						quit()
				if gc_winsize != file_gc_winsize:
					print 'ERROR: GC Content Window Sizes of files to not match!'
					quit()
		elif line.startswith('#GC Content of'):
			gc_loc_counts[gc_linecount] += int(line.split()[-2])
			gc_linecount += 1
		elif not line.startswith('#'):
			infile.close()
			break

outfile = open(args.output, 'w')

outfile.write('#Genome Window Size: %i bp\tGC Content Window Size: %.*f\n' %(genome_winsize, -gc_places, gc_winsize))
for idx, window in enumerate(gc_loc_counts):
	outfile.write('#GC Content of %.*f-%.*f: %i locations\n' %(-gc_places + 2, idx * gc_winsize, -gc_places + 2, min(1., (idx + 1) * gc_winsize - decimal.Decimal((0, (0,1), gc_places - 2))), window))


for infilename in args.infiles:
	infile = open(infilename, 'r')
	for line in infile:
		if not line.startswith('#'):
			outfile.write(line)
	infile.close()

outfile.close()
