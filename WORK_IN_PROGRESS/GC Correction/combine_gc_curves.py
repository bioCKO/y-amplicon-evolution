#!/usr/bin/python

'''
Combines two or more gc curve files.
'''

import os
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio.Statistics import lowess

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Combines two or more gc curve files")
	parser.add_argument("curves", help = "GC curves to combine", nargs="+")
	parser.add_argument("-o", "--output", help = "Output filename of combined GC curve (default: input filenames + .gc_curve)")
	parser.add_argument("-p", "--plot", help = "Plot the combines GC curve and save plot as output filename + .png", action="store_true")

	args = parser.parse_args()
	
	if args.output:
		outname = args.output
	else:
		outname = '_'.join([x.split('/')[-1][:-9] if x.endswith('.gc_curve') 
							else x for x in args.curves]) + '.gc_curve'
		
	if os.path.isfile(outname):
		print 'Output file %s already exists!' %(outname)
		quit()
	
	read_counts = []
	loc_counts = []

	for gcfile in args.curves:
		ind_read_counts = []
		ind_loc_counts = []
		infile = open(gcfile, 'r')
		for line in infile:
			if line.startswith('#'):
				continue
			data = line.split()
			ind_read_counts.append(int(data[3]))
			ind_loc_counts.append(int(data[4]))
		infile.close()
		read_counts.append(np.array(ind_read_counts))
		loc_counts.append(np.array(ind_loc_counts))
		
	read_counts = sum(read_counts)
	loc_counts = sum(loc_counts)
	
	#Normalize to mean coverage in the entire file
	read_norm, loc_norm = sum(read_counts), sum(loc_counts)#(alignment_bam.mapped + alignment_bam.unmapped)/2, sum(alignment_bam.lengths)
	#Have to normalize the old way if I'm not going into the source bam files, which I might add later?
	gc_curve = [0.0 if loc_count == 0 else float(read_count * loc_norm)/(loc_count * read_norm) 
				for read_count, loc_count in zip(read_counts, loc_counts)]
	
	
	#Correct any extreme outliers caused by low read count
	outliers = []
	for window in xrange(len(gc_curve)):
		if read_counts[window] < 10:
			if window == 0 and gc_curve[window] - 0.5 > gc_curve[window + 1]:
				outliers.append(window)
				gc_curve[window] = gc_curve[window + 1]
			elif window == len(gc_curve) - 1 and gc_curve[window] - 0.5 > gc_curve[window - 1]:
				outliers.append(window)
				gc_curve[window] = gc_curve[window + 1]
			elif gc_curve[window] - 0.5 > gc_curve[window - 1] and gc_curve[window] - 0.5 > gc_curve[window + 1]:
				outliers.append(window)
				gc_curve[window] = (gc_curve[window - 1] + gc_curve[window + 1])/2.
	
	gc_x = np.array([x / float(len(gc_curve) - 1) for x in xrange(len(gc_curve))])
	smoothed_gc_curve = lowess.lowess(gc_x, gc_curve, f = 0.1, iter = 1)
	smoothed_gc_curve = [max(0.0, x) for x in smoothed_gc_curve]
	
	outfile = open(outname, 'w')

	outfile.write('# GC Curve file combined from %s\n' %(', '.join(args.curves)))
	outfile.write('# Curve calculated from %i reads at %i locations\n' %(sum(read_counts), sum(loc_counts)))
	if len(outliers) > 0:
		outfile.write('# Windows with corrected extreme outlier GC bias: %s\n' %(', '.join(['%.*f-%.*f' %(len(str((len(gc_curve) - 1))) + 2, window * 1.0/(len(gc_curve) - 1), len(str((len(gc_curve) - 1))) + 2, min(1., (window + 1) * 1.0/(len(gc_curve) - 1) - (1. / 10 ** (len(str((len(gc_curve) - 1))) + 2)))) for window in outliers])))
	outfile.write('#\n')
	outfile.write('#GC_content\tSmoothed_GC_bias\tRaw_GC_Bias\tNo_of_reads\tNo_of_locations\n')
	for window, bias in enumerate(gc_curve):
		outfile.write('%.*f-%.*f\t%f\t%f\t%i\t%i\n' %(len(str((len(gc_curve) - 1))) + 2, window * 1.0/(len(gc_curve) - 1), len(str((len(gc_curve) - 1))) + 2, min(1., (window + 1) * 1.0/(len(gc_curve) - 1) - (1. / 10 ** (len(str((len(gc_curve) - 1))) + 2))), smoothed_gc_curve[window], bias, read_counts[window], loc_counts[window]))

	if args.plot:
		plt.clf()
		plt.plot([x/float((len(gc_curve) - 1)) for x in range((len(gc_curve) - 1) + 1)], gc_curve, 'o')
		plt.plot([x/float((len(gc_curve) - 1)) for x in range((len(gc_curve) - 1) + 1)], smoothed_gc_curve, 'red')
		plt.xlabel('GC content')
		plt.ylabel('GC bias')
		plt.title('GC Bias in %s' %(outname))
		plt.savefig(outname + '.png')

	outfile.close()
	
