#!/usr/bin/python

import pickle
import sys
import os
from make_chrY_vector import bed_to_depthvec

if os.path.isfile(sys.argv[2]):
	print 'Outfile already exists!'
	quit()

rfile = open('/lab/Page_lab-users/lsteitz/Y_repeats.txt', 'r')

y_reps = {}
for line in rfile:
	if line.startswith('#'):
		continue
	data = line.rstrip().split('\t')
	if not data[5] in y_reps:
		y_reps[data[5]] = []
	y_reps[data[5]].append(data)

rfile.close()

source_bed = sys.argv[1] #'chrY_to_all_100bp_cut14.bedgraph'

'''
			   'chrY_to_all_100bp_cut11.bedgraph', 
			   'chrY_to_all_100bp_cut0.bedgraph', 
			   'chrY_to_all_50bp_cut14.bedgraph', 
			   'chrY_to_all_50bp_cut11.bedgraph', 
			   'chrY_to_all_50bp_cut0.bedgraph', 
			   'chrY_to_all_35bp_cut14.bedgraph', 
			   'chrY_to_all_35bp_cut11.bedgraph', 
			   'chrY_to_all_35bp_cut0.bedgraph'
'''


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


outfile = open(sys.argv[2], 'w')
pickle.dump(goodbase_set, outfile)
outfile.close()

