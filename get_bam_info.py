#!/usr/bin/python

from Bio import SeqIO
import pysam
import numpy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_bam_libs(pybamfile):
	libs = []
	for read in pybamfile.fetch():
		if read.qname.split('.')[0] not in libs:
			libs.append(read.qname.split('.')[0])
	return libs

def get_bam_info(pybamfile, reads_per_chr=1000, min_mapq=30):
	#Takes pysam Samfile object as input; returns read length, fragment size, and fragment stddev
	inlens = [] #list of fragment lengths
	readlens = [] #list of read lengths

	for refname in pybamfile.references:
		if refname == 'chrM': #Skip mitochondrial sequence
			continue
		counter = 0
		for read in pybamfile.fetch(refname):
			if read.is_proper_pair and read.is_read1 and read.mapq >= min_mapq and read.alen == read.rlen:
				mate = pybamfile.mate(read)
				if  mate.mapq >= min_mapq and mate.alen == mate.rlen and read.tid == mate.tid: 
				#Use only paired reads where both reads in the pair have mapping quality >= a cutoff (default=30) and both reads fully map to the reference (no soft clipping) on the same chromosome
					inlens.append(abs(read.tlen))
					readlens.append(read.rlen)
					counter += 1
					if counter >= reads_per_chr:
						break
	#Remove top and bottom 2% of fragment lengths, to get rid of extreme outliers that will have a major effect on mean or stddev
	inlens.sort()
	inlens = inlens[int(len(inlens)*.02):-int(len(inlens)*.02)]
	return inlens, readlens

def get_bam_info_per_library(pybamfile, reads_per_chr=1000, liblist=None, min_mapq=30):
	#Takes pysam Samfile object as input; returns read length, fragment size, and fragment stddev
	inlens = {} #list of fragment lengths
	readlens = {} #list of read lengths
	counters = {}
	if liblist:
		for lib in liblist:
			inlens[lib] = []
			readlens[lib] = []
			counters[lib] = 0
	for refname in pybamfile.references:
		if refname == 'chrM': #Skip mitochondrial sequence
			continue
		#print refname
		for lib in counters.keys():
			counters[lib] = 0
		for read in pybamfile.fetch(refname):
			readrun = read.qname.split('.')[0]
			if readrun not in inlens:
				inlens[readrun] = []
				readlens[readrun] = []
				counters[readrun] = 0
			else:
				if counters[readrun] >= reads_per_chr:
					continue
			if read.is_proper_pair and read.is_read1 and read.mapq >= min_mapq and read.alen == read.rlen:
				mate = pybamfile.mate(read)
				if  mate.mapq >= min_mapq and mate.alen == mate.rlen and read.tid == mate.tid: 
				#Use only paired reads where both reads in the pair have mapping quality >= a cutoff (default=30) and both reads fully map to the reference (no soft clipping) on the same chromosome
					inlens[readrun].append(abs(read.tlen))
					readlens[readrun].append(read.rlen)
					counters[readrun] += 1
					if all(count >= reads_per_chr for count in counters.values()):
						break
	#Remove top and bottom 2% of fragment lengths, to get rid of extreme outliers that will have a major effect on mean or stddev
	for library in inlens:
		inlens[library].sort()
		inlens[library] = inlens[library][int(len(inlens[library])*.02):-int(len(inlens[library])*.02)]
	return inlens, readlens

def get_simple_bam_info(pybamfile, no_reads=10000, min_mapq=30):
	inlens = [] #list of fragment lengths
	readlens = [] #list of read lengths
	counter = 0
	for read in pybamfile.fetch():
		if read.is_proper_pair and read.is_read1 and read.mapq >= min_mapq and read.alen == read.rlen:
			mate = pybamfile.mate(read)
			if  mate.mapq >= min_mapq and mate.alen == mate.rlen and read.tid == mate.tid: 
			#Use only paired reads where both reads in the pair have mapping quality >= a cutoff (default=30) and both reads fully map to the reference (no soft clipping) on the same chromosome
				inlens.append(abs(read.tlen))
				readlens.append(read.rlen)
				counter += 1
				if counter >= no_reads:
					break
	#Remove top and bottom 2% of fragment lengths, to get rid of extreme outliers that will have a major effect on mean or stddev
	inlens.sort()
	inlens = inlens[int(len(inlens)*.02):-max(int(len(inlens)*.02),1)]
	return inlens, readlens

def get_frag_hist(in_dist, outname=False):
	fragsizes = numpy.histogram(in_dist, bins=range(min(in_dist)-10, max(in_dist)+10,3))#max(1,(max(in_dist)-min(in_dist))/200))) #The last argument in range can be set constant or to the commented out term, which tries to fit it best to the range of the data
	if outname:
		outtext = open('%s.txt' %(outname), 'w')
		outtext.write('#Fragsize_start\tFragments\n')
		for size, frags in zip(fragsizes[1], fragsizes[0]):
			outtext.write('%i\t%i\n' %(size, frags))
		outtext.close()
	else:
		return fragsizes
	
def plot_frag_size(in_dist, outname, bamfile):
	plt.clf()
	plt.hist(in_dist, bins=range(min(in_dist)-10, max(in_dist)+10,3))#max(1,(max(in_dist)-min(in_dist))/200))) #The last argument in range can be set constant or to the commented out term, which tries to fit it best to the range of the data
	plt.xlabel('Fragment size')
	plt.ylabel('Number of fragments')
	plt.title('Fragment size distribution of %s' %(bamfile))
	plt.savefig('%s.png' %(outname), dpi=250)
	plt.clf()

def calculate_bam_stats(in_dist, read_dist):
	read_len = int(stats.mode(read_dist)[0])
	frag_size = numpy.mean(in_dist)
	frag_dev = numpy.std(in_dist)
	return read_len, frag_size, frag_dev


if __name__ == '__main__':
	
	import argparse
	
	parser = argparse.ArgumentParser(description = "Calculates read and fragment length statistics of a BAM file")
	parser.add_argument("input", help = "Input BAM file used to calculate statistics")
	parser.add_argument("-s", "--simple", help = "Use the first reads, rather than reads from each chromosome, to calculate statistics (overrides --by_library)", action = "store_true")
	parser.add_argument("-l", "--by_library", help = "Return statistics for each library in the BAM file", action = "store_true")
	parser.add_argument("-a", "--all_libraries", help = "Ensures that all libraries are found when --by_library is used, but adds significant amount of time to run", action = "store_true")
	parser.add_argument("-r", "--reads", help = "Number of reads used per chromosome to calculate statistics (or total number of reads with --simple)", type = int)
	parser.add_argument("-p", "--plot", help = "Plot the insert size distribution", action = "store_true")
	parser.add_argument("-o", "--output", help = "Output filename of the plot (default: input filename + _fragsize.png)", default = None)
	parser.add_argument("-q", "--min_mapq", help = "Minimum mapping quality score at which to count pairs (default = 30)", default = 30)
	
	args = parser.parse_args()
	
	bamfile = pysam.Samfile(args.input, 'rb')
	bamname = args.input.split('/')[-1]
	if args.output:
		outname = args.output
	else:
		if args.input.endswith('.bam'):
			outname = '%s_fragsize' %('.'.join(args.input.split('/')[-1].split('.')[:-1]))
		else:
			outname = '%s_fragsize' %(args.input.split('/')[-1])
	
	
	if args.simple:
		if args.reads:
			in_dist, read_dist = get_simple_bam_info(bamfile, no_reads=args.reads, min_mapq=args.min_mapq)
		else:
			in_dist, read_dist = get_simple_bam_info(bamfile, min_mapq=args.min_mapq)
		read_len, frag_size, frag_dev = calculate_bam_stats(in_dist, read_dist)
		print 'Read length: %i\tMean fragment size: %f\tFragment standard deviation: %f' %(read_len, frag_size, frag_dev)	
		get_frag_hist(in_dist, outname=outname)
		if args.plot:
			plot_frag_size(in_dist, outname, bamname)

	
	elif args.by_library:
		if args.all_libraries:
			liblist = get_bam_libs(bamfile)
			if args.reads:
				in_dist, read_dist = get_bam_info_per_library(bamfile, reads_per_chr=args.reads, liblist=liblist, min_mapq=args.min_mapq)
			else:
				in_dist, read_dist = get_bam_info_per_library(bamfile, liblist=liblist, min_mapq=args.min_mapq)
		
		else:
			if args.reads:
				in_dist, read_dist = get_bam_info_per_library(bamfile, reads_per_chr=args.reads, min_mapq=args.min_mapq)
			else:
				in_dist, read_dist = get_bam_info_per_library(bamfile, min_mapq=args.min_mapq)
	
		for library in in_dist:
			read_len, frag_size, frag_dev = calculate_bam_stats(in_dist[library], read_dist[library])
			print 'Library %s\nRead length: %i\tMean fragment size: %f\tFragment standard deviation: %f\n' %(library, read_len, frag_size, frag_dev)
			get_frag_hist(in_dist[library], outname='%s_%s' %(outname, library))
			if args.plot:
				plot_frag_size(in_dist[library], '%s_%s' %(outname, library), '%s (library %s)' %(bamname, library))

		
	else:
		if args.reads:
			in_dist, read_dist = get_bam_info(bamfile, reads_per_chr=args.reads, min_mapq=args.min_mapq)
		else:
			in_dist, read_dist = get_bam_info(bamfile, min_mapq=args.min_mapq)
		read_len, frag_size, frag_dev = calculate_bam_stats(in_dist, read_dist)
		print 'Read length: %i\tMean fragment size: %f\tFragment standard deviation: %f' %(read_len, frag_size, frag_dev)
		get_frag_hist(in_dist, outname=outname)
		if args.plot:
			plot_frag_size(in_dist, outname, bamname)
	
	
	
	
	
	
	
	
	
	
	
	
	
	