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
from scipy.stats import mannwhitneyu
from make_chrY_vector import get_reg_depth


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

def get_reg_wins(depthvec,rstart,rend, mask_repeats=None, winsize=100, masked_cutoff=0.25, skip_masked=True):
	
	wins = []
	for winstart in xrange(0, rend - rstart, winsize):
		if mask_repeats:
			masked_bases = 0.
			for base in xrange(rstart + winstart, min(rstart + winstart + winsize, rend)):
				if not base in mask_repeats:
					masked_bases += 1.
			if masked_bases / winsize >= masked_cutoff:
				if not skip_masked:
					wins.append(0)
				continue
		wins.append(np.mean(depthvec[rstart + winstart:min(rstart + winstart + winsize, rend)]))
	return wins



def get_amplicon_wins(depthvec, rep_dict, *repregs, **kwargs):
	if not 'mask_repeats' in kwargs:
		mask_repeats = None
	else:
		mask_repeats = kwargs['mask_repeats']
		
	if not 'skip_masked' in kwargs:
		skip_masked = True
	else:
		skip_masked = kwargs['skip_masked']
	
	if not 'winsize' in kwargs:
		winsize = 100
	else:
		winsize = kwargs['winsize']
		
	for kwarg in kwargs:
		if kwarg not in ['mask_repeats', 'skip_masked', 'winsize']:
			print 'WARNING: unrecognized keyword argument'
	
	amp_wins = {}
	
	for rep in repregs:
		for repeat_copy in rep_dict[rep]:
			if not repeat_copy[6] in amp_wins:
				amp_wins[repeat_copy[6]] = []
	
			if repeat_copy[3] == 'F':
				amp_wins[repeat_copy[6]] += get_reg_wins(depthvec, 
														 int(repeat_copy[1]),
														 int(repeat_copy[2]),
														 mask_repeats=mask_repeats, 
														 skip_masked=skip_masked, 
														 winsize=winsize)
			elif repeat_copy[3] == 'R':
				amp_wins[repeat_copy[6]] += get_reg_wins(depthvec, 
														 int(repeat_copy[1]),
														 int(repeat_copy[2]),
														 mask_repeats=mask_repeats, 
														 skip_masked=skip_masked, 
														 winsize=winsize)[::-1]
			else:
				print 'ERROR: Amplicon has no orientation'
	
	return amp_wins



if __name__ == '__main__':
	infiles = sys.argv[1:]
#	outname = sys.argv[-1]
# 	if os.path.isfile(outname):
# 		print 'ERROR: Output file already exists.'
# 		quit()
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
	
	normcov, normbases = get_reg_depth(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs)
	normdepth = normcov/normbases
	
	normwins = get_reg_wins(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs)
	
	maxp = [1., 0., 0., 0.]
	for splitloc in xrange(1, len(normwins)-1):
		pval = mannwhitneyu(normwins[:splitloc], normwins[splitloc:]).pvalue
		if pval < maxp[0]:
			maxp = [pval, splitloc, 
					np.mean(normwins[:splitloc]), 
					np.mean(normwins[splitloc:])]
	
	#if abs(maxp[2] - maxp[3])/normdepth > 0.1:		
	sys.stdout.write('Norm reg\t%.4e\t%i\t%.3f\t%.3f\n' %tuple(maxp))
	
	

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

	amplicons = ['IR3', 'TSPY', 'IR1', 'P8', 'P7', 
				 'P6', 'P5', 'IR5', 'P4', 'DYZ19', 
				 'RBMY', 'IR2', 'Blue', 'Teal', 
				 'Green', 'Red', 'DAZ', 'Gray', 'Yellow']

	
	amplicon_regions = {'IR3':['IR3-1', 'IR3-2', 'IR3-3'], 'TSPY':['TSPY'], 
						'IR1':['IR1-1', 'IR1-2', 'IR1-3'], 'P8':['P8'], 
						'P7':['P7'], 'P6':['P6'], 'P5':['P5'], 
						'IR5':['IR5'], 'P4':['P4'], 'DYZ19':['DYZ19'], 
						'RBMY':['RBMY1-1', 'IR2-RBMY1-spacer', 'RBMY1-2'], 
						'IR2':['IR2'], 'Blue':['Blue', 'Blue-plus'], 
						'Teal':['Teal-1', 'Teal-2'], 'Green':['Green'], 
						'Red':['Red-1', 'Red-2'], 'DAZ':['DAZ'], 'Gray':['Gray'], 
						'Yellow':['Yellow-1', 'P1.1/2', 'P1.1/2-spacer', 
							'P1-chr15-1', 'P1.3/4', 'P1.3/4-spacer', 'P1-ch15-2']}
	
	
	unmasked_amps = ['TSPY', 'DYZ19']
	
	copy_numbers = {'Blue': 4,
					'Gray': 2,
					'Green': 3,
					'IR1': 2,
					'IR2': 2,
					'IR3': 2,
					'IR5': 4,
					'P4': 2,
					'P5': 2,
					'P6': 2,
					'P7': 2,
					'P8': 2,
					'Red': 4,
					'Teal': 2,
					'Yellow': 2}
 	
	for amplicon in amplicons:
		if amplicon in unmasked_amps:
			continue
			#amp_wins = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[amplicon])
		elif amplicon not in copy_numbers:
			continue
		else:
			amp_wins = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[amplicon], mask_repeats=goodlocs)
		
		maxps = {}
		
		#This version calculates all split point for all amplicons
		for ampcopy in amp_wins:
			maxp = [1., 0, 0., 0.]
			if len(set(amp_wins[ampcopy])) == 1:
				maxp = [1., 0, amp_wins[ampcopy][0], amp_wins[ampcopy][0]]
				r1cn = int((maxp[2]/normdepth) * copy_numbers[amplicon] + 0.5)
				r2cn = int((maxp[3]/normdepth) * copy_numbers[amplicon] + 0.5)
				maxps[ampcopy] = maxp + [r1cn, r2cn]
				continue
			for splitloc in reversed(xrange(1, len(amp_wins[ampcopy])-1)):
				pval = mannwhitneyu(amp_wins[ampcopy][:splitloc], amp_wins[ampcopy][splitloc:]).pvalue
				if pval < maxp[0]:
					maxp = [pval, splitloc, 
							np.mean(amp_wins[ampcopy][:splitloc]), 
							np.mean(amp_wins[ampcopy][splitloc:])]
			
			r1cn = int((maxp[2]/normdepth) * copy_numbers[amplicon] + 0.5)
			r2cn = int((maxp[3]/normdepth) * copy_numbers[amplicon] + 0.5)
			maxps[ampcopy] = maxp + [r1cn, r2cn]			
			
		sys.stdout.write('%s' %(amplicon))
		for ampreg in maxps:
			sys.stdout.write('\t%s\t%.4e\t%i\t%.3f\t%.3f\t%i\t%i' %tuple([ampreg] + maxps[ampreg]))
						
		sys.stdout.write('\n')
		'''
		#This version stops calculating if an amplicon won't fulfill the criteria for reporting
		for ampcopy in amp_wins:
			maxp = [1., 0, 0., 0.]
			for splitloc in xrange(1, len(amp_wins[ampcopy])-1):
				pval = mannwhitneyu(amp_wins[ampcopy][:splitloc], amp_wins[ampcopy][splitloc:]).pvalue
				if pval < maxp[0]:
					maxp = [pval, splitloc, 
							np.mean(amp_wins[ampcopy][:splitloc]), 
							np.mean(amp_wins[ampcopy][splitloc:])]
			
			r1cn = int((maxp[2]/normdepth) * copy_numbers[amplicon] + 0.5)
			r2cn = int((maxp[3]/normdepth) * copy_numbers[amplicon] + 0.5)
			
			if r1cn == r2cn: #Different CN in all copies
				break
			#if maxp[1] < 25 or maxp[1] > len(amp_wins[ampcopy]) - 25: #Breakpoint isn't too close to an edge
			#	break
			
			maxps[ampcopy] = maxp + [r1cn, r2cn]
			if len(maxps) == len(amp_wins):
				if 1: #max([x[1] for x in maxps.values()]) - min([x[1] for x in maxps.values()]) <= 25: #Breakpoint is in similar location in each copy
					sys.stdout.write('%s' %(amplicon))
					for ampreg in maxps:
						sys.stdout.write('\t%s\t%.4e\t%i\t%.3f\t%.3f\t%i\t%i' %tuple([ampreg] + maxps[ampreg]))
						
					sys.stdout.write('\n')
					sys.stdout.flush()
		'''
	
	#Control regions
	ctrl_regs = [[4000000, 5000000], 
				 [7000000, 7500000], 
				 [8000000, 8500000], 
				 [16780000, 16800000], 
				 [19550000, 19560000], 
				 [20500000, 20550000]]
	
	for cridx, ctrl_reg in enumerate(ctrl_regs):
		ctrl_wins = get_reg_wins(chrY_depth, ctrl_reg[0], ctrl_reg[1], mask_repeats=goodlocs)
		maxp = [1., 0., 0., 0.]
		for splitloc in xrange(1, len(ctrl_wins)-1):
			pval = mannwhitneyu(ctrl_wins[:splitloc], ctrl_wins[splitloc:]).pvalue
			if pval < maxp[0]:
				maxp = [pval, splitloc, 
						np.mean(ctrl_wins[:splitloc]), 
						np.mean(ctrl_wins[splitloc:])]
		
		#if abs(maxp[2] - maxp[3])/normdepth > 0.1:
		sys.stdout.write('Ctrl_reg_%i\t%.4e\t%i\t%.3f\t%.3f\n' %tuple([cridx] + maxp))
		
	sys.stdout.flush()		

	'''
	
	regnames = ['IR3', 'TSPY', 'IR1', 'P8', 'P7', 'P6', 'P5', 'IR5', 'P4', 'DYZ19', 'RBMY', 'IR2', 'Blue', 'Teal', 'Green', 'Red', 'DAZ', 'Gray', 'Yellow']
	regnames += ['Ctrl_reg_%i' %(x) for x in xrange(len(ctrl_regs))]
	
	meds = [IR3_med, TSPY_med, IR1_med, P8_med, P7_med, P6_med, P5_med, IR5_med, P4_med, DYZ19_med,
			RBMY_med, IR2_med, Blue_med, Teal_med, Green_med, Red_med, DAZ_med, Gray_med, Yellow_med]
	meds += ctrl_meds
	
	normed_meds = [x/norm_med for x in meds]
	'''