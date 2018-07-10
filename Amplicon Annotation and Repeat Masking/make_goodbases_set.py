#!/usr/bin/python

'''
Takes bedGraph file of Y-chromosome mapability created by make_rep_graph.py and 
creates a python pickle file of the set of bases that map as many times as expected in 
the amplicons and control regions used in Teitz et al. 2018.
'''

import pickle
import sys
import os
from make_chrY_vector import bed_to_depthvec

source_bed = sys.argv[1] #bedgraph from make_rep_graph.py
outfile = sys.argv[2] #name of output pickle file

if os.path.isfile(outfile):
	print 'Outfile already exists!'
	quit()

rfile = open('./Y_repeats.txt', 'r')

y_reps = {}
for line in rfile:
	if line.startswith('#'):
		continue
	data = line.rstrip().split('\t')
	if not data[5] in y_reps:
		y_reps[data[5]] = []
	y_reps[data[5]].append(data)

rfile.close()


y_rep_pcts = {x:[] for x in y_reps}
winlen = int(source_bed.split('/')[-1].split('_')[3][:-2])
bases = bed_to_depthvec(source_bed, 'chrY')

goodbase_set = set([])

for reg in y_reps:
	for subreg in y_reps[reg]:	
		for base in xrange(int(subreg[1]),int(subreg[2])):
			if all([0 < x <= int(subreg[8]) for x in bases[base - winlen + 1:base + 1]]):
				goodbase_set.add(base)

#Control regions
controls = [[14500000, 15500000], [4000000, 5000000], 
				 [7000000, 7500000], 
				 [8000000, 8500000], 
				 [16780000, 16800000], 
				 [19550000, 19560000], 
				 [20500000, 20550000]]

for control in controls:
	for base in xrange(control[0], control[1]):
		if all([x == 1 for x in bases[base - winlen + 1:base + 1]]):
			goodbase_set.add(base)


outfile = open(outfile, 'w')
pickle.dump(goodbase_set, outfile)
outfile.close()

