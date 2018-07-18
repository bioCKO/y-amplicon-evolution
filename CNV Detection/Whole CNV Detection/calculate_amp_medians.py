#!/usr/bin/python

'''
Takes bedgraph coverage file(s) and returns a file of median normalized depths of chrY amplicons and control regions
'''

import os
import argparse
import pickle
import numpy as np

def bed_to_depthvec(bedfile, chrom, regfile=None):
	depthvec = []
	last = -1
	with open(bedfile, 'r') as infile:
		for line in infile:
			data = line.split()
			if data[0] == chrom:
				depthvec += [0] * (int(data[1]) - last - 1)
				for thing in xrange(int(data[1]), int(data[2])):
					depthvec.append(float(data[3]))
					last = thing
	if regfile:
		with open(regfile, 'r') as infile:
			goodregs = []
			for line in infile:
				data = line.split()
				goodregs.append([int(data[1]), int(data[2])])
		tempvec = []
		for reg in goodregs:
			tempvec += depthvec[reg[0]:reg[1]]
		depthvec = tempvec
	return depthvec

def get_reg_wins(depthvec,rstart,rend, mask_repeats=None, winsize=100):
	depths = []
	for base, depth in enumerate(depthvec[rstart:rend]):
		if mask_repeats and base + rstart in mask_repeats:
			depths.append(depth)
		elif not mask_repeats:
			depths.append(depth)
	wins = []
	for winstart in xrange(0, len(depths), winsize):
		wins.append(sum(depths[winstart:winstart + winsize])/float(len(depths[winstart:winstart + winsize])))
	return wins
	

def get_repeat_wins(depthvec, rep_dict, rname, buffer=10, mask_repeats=None):
	wins = []
	for repeat_copy in rep_dict[rname]:
		wins += get_reg_wins(depthvec, 
										int(repeat_copy[1]) + buffer,
										int(repeat_copy[2]) - buffer,
										mask_repeats=mask_repeats)
	return wins

def get_multi_rep_wins(depthvec, rep_dict, *repregs, **kwargs):
	if not 'mask_repeats' in kwargs:
		mask_repeats = None
	else:
		mask_repeats = kwargs['mask_repeats']
	if 'buffer' in kwargs:
		buffer = kwargs['buffer']
	else:
		buffer = 10
	for kwarg in kwargs:
		if kwarg not in ['mask_repeats', 'buffer']:
			print 'WARNING: unrecognized keyword argument'
	wins = []
	for rep in repregs:
		wins += get_repeat_wins(depthvec, rep_dict, rep, buffer=buffer, mask_repeats=mask_repeats)
	return wins


if __name__ == '__main__':


	dirname = '/'.join(os.path.realpath(__file__).split('/')[:-3])

	parser = argparse.ArgumentParser(description = "Takes bedgraph coverage file(s) and returns a file of median normalized depths of chrY amplicons and control regions")
	
	parser.add_argument("input_files", nargs = '*', help = "bedGraph file(s) of Y chromosome depth")
	parser.add_argument("-o", "--outfile", default = "./amplicon_medians", help = "Output filename")
	parser.add_argument("-g", "--good_locations", 
						default = "%s/Amplicon Annotation and Repeat Masking/good_chrY_amplicon_bases_w100_c11.pickle" %(dirname), 
						help = "File of unmasked Y chromosome locations (provided in Amplicon Annotation and Repeat Masking")
	parser.add_argument("-r", "--repeats_file", 
						default = "%s/Amplicon Annotation and Repeat Masking/Y_repeats.txt" %(dirname), 
						help = "File of Y chromosome amplicon locations (provided in Amplicon Annotation and Repeat Masking")
	
	args = parser.parse_args()


	if os.path.isfile(args.outfile):
		print 'ERROR: Output file already exists.'
		quit()
	goodlocsfile = open(args.good_locations, 'r')
	goodlocs = pickle.load(goodlocsfile)
	goodlocsfile.close()
	
	chrY_depth = []
	for infile in args.input_files:
		chrY_depth.append(bed_to_depthvec(infile, 'chrY'))
	
	if len(chrY_depth) == 1:
		chrY_depth = chrY_depth[0]
	else:
		max_len = max([len(x) for x in chrY_depth])
		for depth_idx in xrange(len(chrY_depth)):
			chrY_depth[depth_idx] += [0] * (max_len - len(chrY_depth[depth_idx]))
		chrY_depth = [sum(x) for x in zip(*chrY_depth)]
	
	#Get median depth at single-copy region for normalization
	
	normwins = get_reg_wins(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs)
	norm_med = np.median(normwins)
	
	#Get median coverages at multimapping regions

	rfile = open(args.repeats_file, 'r')

	y_reps = {}
	for line in rfile:
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if not data[5] in y_reps:
			y_reps[data[5]] = []
		y_reps[data[5]].append(data)

	rfile.close()
	
	
	amplicons = ['IR3', 'TSPY', 'IR1', 'P8', 'P7', 
				 'P6', 'P5', 'IR5', 'P4', 'DYZ19', 
				 'RBMY', 'IR2', 'Blue', 'Teal', 
				 'Green', 'Red', 'DAZ', 'Gray', 'Yellow']

	
	amplicon_regions = {'IR3':['IR3-1', 'IR3-2', 'IR3-3'], 'TSPY':['TSPY'], 
						'IR1':['IR1-1', 'IR1-2', 'IR1-3'], 'P8':['P8'], 
						'P7':['P7'], 'P6':['P6'], 'P5':['P5'], 
						'IR5':['IR5'], 'P4':['P4'], 'DYZ19':['DYZ19'], 
						'RBMY':['RBMY1-1', 'IR2-RBMY1-spacer', 'RBMY1-2'], 
						'IR2':['IR2'], 'Blue':['Blue'], 
						'Teal':['Teal-1', 'Teal-2'], 'Green':['Green'], 
						'Red':['Red-1', 'Red-2'], 'DAZ':['DAZ'], 'Gray':['Gray'], 
						'Yellow':['Yellow-1', 'P1.1/2', 'P1.1/2-spacer', 
							'P1-chr15-1', 'P1.3/4', 'P1.3/4-spacer', 'P1-ch15-2']}
	
	
	unmasked_amps = ['TSPY', 'DYZ19']

	medians = []
	
	for amplicon in amplicons:
		if amplicon in unmasked_amps:
			amp_wins = get_multi_rep_wins(chrY_depth, y_reps, *amplicon_regions[amplicon])
		else:
			amp_wins = get_multi_rep_wins(chrY_depth, y_reps, *amplicon_regions[amplicon], mask_repeats=goodlocs)
		medians.append(np.median(amp_wins))


	#Control regions
	ctrl_meds = []
	ctrl_regs = [[4000000, 5000000], 
				 [7000000, 7500000], 
				 [8000000, 8500000], 
				 [16780000, 16800000], 
				 [19550000, 19560000], 
				 [20500000, 20550000]]
	
	for ctrl_reg in ctrl_regs:
		ctrl_wins = get_reg_wins(chrY_depth, ctrl_reg[0], ctrl_reg[1], mask_repeats=goodlocs)
		medians.append(np.median(ctrl_wins))
		
	regnames = amplicons + ['Ctrl_reg_%i' %(x) for x in xrange(len(ctrl_regs))]
	
	normed_meds = [x/norm_med for x in medians]

	outfile = open(args.outfile, 'w')
	outfile.write('#Source file: %s\tSingle-copy median depth: %f\n' %(', '.join(args.input_files), norm_med))
	outfile.write('#%s\n' %('\t'.join(regnames)))
	outfile.write('%s\n' %('\t'.join([str(x) for x in normed_meds])))
	outfile.close()
	