#!/usr/bin/python

'''
Take bedgraph coverage file and return vector of depths of chrY multimapping regions
'''

import pickle
import sys
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def get_reg_wins(depthvec,rstart,rend, mask_repeats=None, winsize=30):
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
	infiles = sys.argv[1:-1]
	outname = sys.argv[-1]
	if os.path.isfile(outname):
		print 'ERROR: Output file already exists.'
		quit()
	#goodlocsfile = open('/lab/solexa_page/lsteitz/chrY_repeatmasker_set.pickle', 'r')
	goodlocsfile = open('/lab/solexa_page/lsteitz/Y_repetition/good_chrY_amplicon_bases_w100_c11.pickle', 'r')
	goodlocs = pickle.load(goodlocsfile)
	goodlocsfile.close()
	
	chrY_depth = []
	for infile in infiles:
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
	plt.clf()
	plt.hist(normwins, bins=30)
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s_normreg.png' %(outname))
	
	
	#norm_med = np.median(normwins)
	
	#Get median coverages at multimapping regions

	rfile = open('/lab/Page_lab-users/lsteitz/Y_repeats.txt', 'r')

	y_reps = {}
	for line in rfile:
		if line.startswith('#'):
			continue
		data = line.split('\t')
		if not data[5] in y_reps:
			y_reps[data[5]] = []
		y_reps[data[5]].append(data)

	rfile.close()

	'''
	#IR3
	IR3_wins = get_multi_rep_wins(chrY_depth, y_reps, 'IR3-1', 'IR3-2', 'IR3-3', mask_repeats=goodlocs)
	IR3_med = np.median(IR3_wins)
	
	#TSPY
	TSPY_wins = get_multi_rep_wins(chrY_depth, y_reps, 'TSPY')
	TSPY_med = np.median(TSPY_wins)

	#IR1
	IR1_wins = get_multi_rep_wins(chrY_depth, y_reps, 'IR1-1', 'IR1-2', 'IR1-3', mask_repeats=goodlocs)
	IR1_med = np.median(IR1_wins)
	
	#P8
	P8_wins = get_multi_rep_wins(chrY_depth, y_reps, 'P8', mask_repeats=goodlocs)
	P8_med = np.median(P8_wins)
	
	#P7
	P7_wins = get_multi_rep_wins(chrY_depth, y_reps, 'P7', mask_repeats=goodlocs)
	P7_med = np.median(P7_wins)

	#P6
	P6_wins = get_multi_rep_wins(chrY_depth, y_reps, 'P6', mask_repeats=goodlocs)
	P6_med = np.median(P6_wins)
	
	#P5
	P5_wins = get_multi_rep_wins(chrY_depth, y_reps, 'P5', mask_repeats=goodlocs)
	P5_med = np.median(P5_wins)

	#IR5
	IR5_wins = get_multi_rep_wins(chrY_depth, y_reps, 'IR5', mask_repeats=goodlocs)
	IR5_med = np.median(IR5_wins)

	#P4
	P4_wins = get_multi_rep_wins(chrY_depth, y_reps, 'P4', mask_repeats=goodlocs)
	P4_med = np.median(P4_wins)

	#DYZ19
	DYZ19_wins = get_multi_rep_wins(chrY_depth, y_reps, 'DYZ19')
	DYZ19_med = np.median(DYZ19_wins)

	#RBMY
	RBMY_wins = get_multi_rep_wins(chrY_depth, y_reps, 'RBMY1-1', 'IR2-RBMY1-spacer', 'RBMY1-2', mask_repeats=goodlocs)
	RBMY_med = np.median(RBMY_wins)

	#IR2
	IR2_wins = get_multi_rep_wins(chrY_depth, y_reps, 'IR2', mask_repeats=goodlocs)
	IR2_med = np.median(IR2_wins)

	#Blue
	Blue_wins = get_multi_rep_wins(chrY_depth, y_reps, 'Blue', 'Blue-plus', mask_repeats=goodlocs)
	Blue_med = np.median(Blue_wins)

	#Teal
	Teal_wins = get_multi_rep_wins(chrY_depth, y_reps, 'Teal-1', 'Teal-2', mask_repeats=goodlocs)
	Teal_med = np.median(Teal_wins)
	'''
	#Green
	Green_wins = get_multi_rep_wins(chrY_depth, y_reps, 'Green', mask_repeats=goodlocs)
	#Green_med = np.median(Green_wins)
	plt.clf()
	plt.hist(Green_wins, bins=30)
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s_green.png' %(outname))
	
	#Red
	Red_wins = get_multi_rep_wins(chrY_depth, y_reps, 'Red-1', 'Red-2', mask_repeats=goodlocs)
	#Red_med = np.median(Red_wins)
	plt.clf()
	plt.hist(Red_wins, bins=30)
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s_red.png' %(outname))
	'''
	#DAZ
	DAZ_wins = get_multi_rep_wins(chrY_depth, y_reps, 'DAZ', mask_repeats=goodlocs)
	DAZ_med = np.median(DAZ_wins)

	#Gray
	Gray_wins = get_multi_rep_wins(chrY_depth, y_reps, 'Gray', mask_repeats=goodlocs)
	Gray_med = np.median(Gray_wins)

	#Yellow
	Yellow_wins = get_multi_rep_wins(chrY_depth, y_reps, 
		'Yellow-1', 'P1.1/2', 'P1.1/2-spacer',
		'P1.1/2', 'P1-chr15-1', 'P1.3/4',
		'P1.3/4-spacer', 'P1.3/4', 'P1-ch15-2', mask_repeats=goodlocs)
	Yellow_med = np.median(Yellow_wins)
	
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
		ctrl_meds.append(np.median(ctrl_wins))
		
	regnames = ['IR3', 'TSPY', 'IR1', 'P8', 'P7', 'P6', 'P5', 'IR5', 'P4', 'DYZ19', 'RBMY', 'IR2', 'Blue', 'Teal', 'Green', 'Red', 'DAZ', 'Gray', 'Yellow']
	regnames += ['Ctrl_reg_%i' %(x) for x in xrange(len(ctrl_regs))]
	
	meds = [IR3_med, TSPY_med, IR1_med, P8_med, P7_med, P6_med, P5_med, IR5_med, P4_med, DYZ19_med,
			RBMY_med, IR2_med, Blue_med, Teal_med, Green_med, Red_med, DAZ_med, Gray_med, Yellow_med]
	meds += ctrl_meds
	
	normed_meds = [x/norm_med for x in meds]

	outfile = open(outname, 'w')
	outfile.write('#Source file: %s\tSingle-copy median depth: %f\n' %(', '.join(infiles), norm_med))
	outfile.write('#%s\n' %('\t'.join(regnames)))
	outfile.write('%s\n' %('\t'.join([str(x) for x in normed_meds])))
	outfile.close()
	'''