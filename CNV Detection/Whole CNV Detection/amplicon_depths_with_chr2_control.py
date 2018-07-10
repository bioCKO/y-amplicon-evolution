#!/usr/bin/python

'''
Takes a default Y chromosome amplicon depth file normalized by single-copy Y sequence 
and outputs a amplicon depth file normalized by the depth of chr2 
(used to verify Y-specificity of amplicons with data from XX and XY individuals)
'''

import sys
from bed_to_vector import bed_to_depthvec
from make_chrY_vector import get_reg_depth
import glob
import os
import numpy as np

chrY_vector = sys.argv[1] #filename of Y chromosome amplicon depth file
chr2_bedgraph = sys.argv[2] #bedGraph of chr2 depth from the same individual
outname = sys.argv[3] #output filename
if os.path.isfile(outname):
	print 'Error: Output file already exists'
	quit()

control_regs = [[80000000, 85000000]]

chr2_depthvec = bed_to_depthvec(chr2_bedgraph, 'chr2')
print len(chr2_depthvec)
ctrl_depths = []
for reg in control_regs:
	reg_depth, reg_len = get_reg_depth(chr2_depthvec, reg[0], reg[1])
	ctrl_depths.append(reg_depth/reg_len)

print ctrl_depths

y_vec = open(chrY_vector, 'r')
sc_depth = float(y_vec.readline().rstrip().split()[-1])
amp_names = y_vec.readline()
amp_depths = [float(x) for x in y_vec.readline().rstrip().split()]
amp_depths = [x * sc_depth / np.mean(ctrl_depths) for x in amp_depths]

outfile = open(outname, 'w')
outfile.write('#Source files: %s, %s\tchr2 control depth: %f\n' %(chrY_vector, chr2_bedgraph, np.mean(ctrl_depths)))
outfile.write(amp_names)
outfile.write('\t'.join(map(str, amp_depths)))
outfile.close()

