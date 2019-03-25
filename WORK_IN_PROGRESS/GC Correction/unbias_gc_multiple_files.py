#!/usr/bin/python

'''
Takes one or more bam files and one gc_curve file as input and returns a single bedgraph file 
of de-GC biased coverage.
'''

def read_cigar(cigarstring):
	numb = ''
	cigar = []
	for char in cigarstring:
		if char.isdigit():
			numb += char
		else:
			cigar.append([int(numb),char])
			numb = ''
	return cigar

def len_cigar(cigarstring):
	'Returns template length that read aligns to using cigar string'
	cigar = read_cigar(cigarstring)
	cig_len = 0
	for char in cigar:
		if char[1] == 'M' or char[1] == 'D':
			cig_len += char[0]
	return cig_len			

def get_chromtext(chromfile):
	fafile = open(chromfile, 'r')
	chromtext = fafile.read()
	chromtext = chromtext.split()[1:]
	chromtext = ''.join(chromtext)
	fafile.close()
	return chromtext

if __name__ == '__main__':
	
	import argparse
	import numpy
	from scipy import stats
	import pysam
	import subprocess
	import get_bam_info
	from calculate_gc_curve import gc_content
	from Bio import SeqIO
	
	
	parser = argparse.ArgumentParser(description = "Takes bam file and gc_curve file as input and returns bedgraph file of de-GC biased coverage")
	parser.add_argument("gcfile", help = "gc_curve file generated by calculate_gc_curve")
	parser.add_argument("reference", help = "FASTA file of reference chromosome the BAM file was aligned to")
	parser.add_argument("bamfiles", help = "BAM files to calculate unbiased depth from", nargs = '+')
	parser.add_argument("-o", "--output", help = "Output filename (default: input BAM file name + .ungcbiased.bedgraph)")
	parser.add_argument("-s", "--stats_info", help = ".stats file containing read and fragment length info", nargs = '*')
	parser.add_argument("-c", "--chrom", help = "Chromosome to return adjusted depth of (default = chrY)", default = 'chrY')
	
	args = parser.parse_args()	
	
	if args.stats_info and len(args.stats_info) != len(args.bamfiles):
		print "ERROR: Please enter the same number of BAM files and stats files."
		quit()
	
	#Step 1: Read GC curve
	
	gc_curve = []
	with open(args.gcfile, 'r') as infile:
		for line in infile:
			if line.startswith('#'):
				continue
			gc_curve.append(float(line.split()[1]))
	gc_bins = len(gc_curve) -1 #The -1 account for the line for 100% GC content
	
	#Step 2: Make empty depth vector
	
	chromref = get_chromtext(args.reference)
	chromvect = [0] * len(chromref)
	
	#Step 3: Get depth of individual BAM files
	
	for input_idx, bamfilename in enumerate(args.bamfiles):
		bamfile = pysam.Samfile(bamfilename, 'rb')

		#Step 3a: Get BAM information
	
		if args.stats_info:
			with open(args.stats_info[input_idx], 'r') as infile:
				bamstats = infile.readline().split()
				read_len = int(bamstats[2])
				frag_size = float(bamstats[6])
				frag_dev = float(bamstats[10])
		else:
			frag_lens, read_lens = get_bam_info.get_bam_info(bamfile)
			read_len = int(stats.mode(read_lens)[0])
			frag_size = numpy.mean(frag_lens)
			frag_dev = numpy.std(frag_lens)
		
		#Step 3b: Build depth vector
		
		alignment = subprocess.Popen('samtools view %s %s' %(bamfilename, args.chrom), shell=True, stdout=subprocess.PIPE, bufsize = 1)
		for readline in alignment.stdout:
			read = readline.split()
			if read[8] == '0' and read[5] == '*':
				continue
			ref_bases = []
			cigar = read_cigar(read[5])
			gpos = 0
			for ciggroup in cigar:
				if ciggroup[1] == 'M':
					ref_bases.append([gpos, gpos + ciggroup[0]])
					gpos += ciggroup[0]
				elif ciggroup[1] == 'D':
					gpos += ciggroup[0]
		
			if read[6] == '=' and 0 < abs(int(read[8])) < frag_size + frag_dev * 10:
				if int(read[8]) > 0:
					gc_cont = gc_content(chromref[int(read[3]) - 1:int(read[3]) + int(read[8]) - 1], int(read[8]))
				elif int(read[8]) < 0:
					gc_cont = gc_content(chromref[int(read[3]) + gpos + int(read[8]) - 1:int(read[3]) + gpos - 1], -int(read[8]))
			else:
				gc_cont = gc_content(read[9], len(read[9]))
		
			readvalue = min(3, 1 / (gc_curve[int(gc_cont * gc_bins)] + .0000001))
				
			for span in ref_bases:
				for pos in xrange(int(read[3]) + span[0] - 1, int(read[3]) + span[1] - 1):
					chromvect[pos] += readvalue
		
		bamfile.close()
	
	#Step 4: Write bedgraph file
	
	if args.output:
		outname = args.output
	else:
		if bamfile.endswith('bam') and '.' in bamfile:
			outname = '%s_%s.ungcbiased.bedgraph' %('.'.join(bamfile.split('.')[:-1]), args.chrom)
		else:
			outname = '%s_%s.ungcbiased.bedgraph' %(bamfile, args.chrom)

	outfile = open(outname, 'w')
	
	startdepth = 0
	startpos = 0
	for pos, depth in enumerate(chromvect):
		if depth != startdepth:
			if startdepth != 0:
				outfile.write('%s\t%i\t%i\t%f\n' %(args.chrom, startpos, pos, startdepth))
			startdepth = depth
			startpos = pos
	
	outfile.close()
	




















	
	