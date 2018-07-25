#!/usr/bin/python

'''
Takes bedgraph coverage file and returns the maximally significant changepoint in depth of each 
amplicon using a modified version of the binary segmentation algorithm. Also plots amplicon 
depth by 100-bp windows and shows location of maximally significant changepoint.
'''

import pickle
import argparse
import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


def get_reg_depth(depthvec,rstart,rend, mask_repeats=None):
	cov = 0
	bases = 0
	for base, depth in enumerate(depthvec[rstart:rend]):
		if mask_repeats and base + rstart in mask_repeats:
			cov += depth
			bases += 1
		elif not mask_repeats:
			cov += depth
			bases += 1
	return cov, float(bases)

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
					wins.append(-1)
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


def get_loc_corrs(winstart, winend, mask_file):
	#Masking bases means that some windows will be skipped; this corrects the window locations that count 
	#only unskipped windows to count all windows
	loc_corrs = {}
	noskipwins = get_reg_wins([0] * 27000000, winstart, winend, mask_repeats=mask_file, skip_masked=False)
	skipdex = 0
	for windex, win in enumerate(noskipwins):
		if win != -1:
			loc_corrs[skipdex] = windex
			skipdex += 1
	return loc_corrs
	
def get_amp_corrs(rep_dict, mask_file, *ampregs):
	#Performs get_loc_corrs for multi_region amplicons
	loc_corrs = {}
	noskipwins = get_amplicon_wins([0] * 27000000, rep_dict, *ampregs, mask_repeats=mask_file, skip_masked=False)
	
	for ampcopy in noskipwins:
		loc_corrs[ampcopy] = {}
		skipdex = 0
		for windex, win in enumerate(noskipwins[ampcopy]):
			if win != -1:
				loc_corrs[ampcopy][skipdex] = windex
				skipdex += 1
	
	return loc_corrs



if __name__ == '__main__':

	dirname = '/'.join(os.path.realpath(__file__).split('/')[:-3])
	
	parser = argparse.ArgumentParser(description = "")
	
	parser.add_argument("input_files", nargs = "*", help = "bedGraph file(s) of Y chromosome depth")
	parser.add_argument("-o", "--outfile", default = "./amplicon_partials", help = "Output filename")
	parser.add_argument("-s", "--segmentation", action = "store_true", help = "Use binary segmentation to find maximum likelihood breakpoint of amplicons (at least one of -s and -p must be chosen)")
	parser.add_argument("-p", "--plot_depths", nargs = "*", help = "create plots of amplicon depth (input is e.g. Blue, Teal, etc.)")
	parser.add_argument("-g", "--good_locations", 
						default = "%s/Amplicon Annotation and Repeat Masking/good_chrY_amplicon_bases_w100_c11.pickle" %(dirname), 
						help = "File of unmasked Y chromosome locations (provided in Amplicon Annotation and Repeat Masking")
	parser.add_argument("-r", "--repeats_file", 
						default = "%s/Amplicon Annotation and Repeat Masking/Y_repeats.txt" %(dirname), 
						help = "File of Y chromosome amplicon locations (provided in Amplicon Annotation and Repeat Masking")

	
	args = parser.parse_args()


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
	
	normcov, normbases = get_reg_depth(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs)
	normdepth = normcov/normbases
	
	
	#Get binary segmentation brakpoints of normalization region:
	if args.segmentation:
		outfile = open(args.outfile, 'w')
		
		normwins = get_reg_wins(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs)
		loc_corrs = get_loc_corrs(14500000, 15500000, goodlocs)
	
		maxp = [1., 0., 0., 0.]
		for splitloc in xrange(1, len(normwins)-1):
			pval = mannwhitneyu(normwins[:splitloc], normwins[splitloc:]).pvalue
			if pval < maxp[0]:
				maxp = [pval, loc_corrs[splitloc], 
						np.mean(normwins[:splitloc]), 
						np.mean(normwins[splitloc:])]
	
		outfile.write('Norm_reg\t%.4e\t%i\t%.3f\t%.3f\n' %tuple(maxp))	

		#Control regions
		ctrl_regs = [[4000000, 5000000], 
					 [7000000, 7500000], 
					 [8000000, 8500000], 
					 [20500000, 20550000]]
	
		for cridx, ctrl_reg in enumerate(ctrl_regs):
			loc_corrs = get_loc_corrs(ctrl_reg[0], ctrl_reg[1], goodlocs)
			ctrl_wins = get_reg_wins(chrY_depth, ctrl_reg[0], ctrl_reg[1], mask_repeats=goodlocs)
			maxp = [1., 0., 0., 0.]
			for splitloc in xrange(1, len(ctrl_wins)-1):
				pval = mannwhitneyu(ctrl_wins[:splitloc], ctrl_wins[splitloc:]).pvalue
				if pval < maxp[0]:
					maxp = [pval, loc_corrs[splitloc], 
							np.mean(ctrl_wins[:splitloc]), 
							np.mean(ctrl_wins[splitloc:])]
		
			outfile.write('Ctrl_reg_%i\t%.4e\t%i\t%.3f\t%.3f\n' %tuple([cridx + 1] + maxp))
	
	
	#Amplicons
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
	
		if (not args.segmentation) and (not amplicon in args.plot_depths):
			continue	
		if amplicon in unmasked_amps:
			continue
		elif amplicon not in copy_numbers:
			continue

		amp_corrs = get_amp_corrs(y_reps, goodlocs, *amplicon_regions[amplicon])
		amp_wins = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[amplicon], mask_repeats=goodlocs)
		
		
		if args.segmentation:
			maxps = {}
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
						maxp = [pval, amp_corrs[ampcopy][splitloc], 
								np.mean(amp_wins[ampcopy][:splitloc]), 
								np.mean(amp_wins[ampcopy][splitloc:])]
			
				r1cn = int((maxp[2]/normdepth) * copy_numbers[amplicon] + 0.5)
				r2cn = int((maxp[3]/normdepth) * copy_numbers[amplicon] + 0.5)
				maxps[ampcopy] = maxp + [r1cn, r2cn]			
			
			outfile.write('%s' %(amplicon))
			for ampcopy in maxps:
				outfile.write('\t%s\t%.4e\t%i\t%.3f\t%.3f\t%i\t%i' %tuple([ampcopy] + maxps[ampcopy]))
						
			outfile.write('\n')
		
		if amplicon in args.plot_depths:
			
			plt.figure(figsize=[7.5, 1.5 * copy_numbers[amplicon]])
			amp_wins_with_masked = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[amplicon], mask_repeats=goodlocs, skip_masked=False)
			max_yval = int(max([np.mean(amp_wins[x]) for x in amp_wins]) + 3 * max([np.std(amp_wins[x]) for x in amp_wins]) + 1)
			
			for ampcopy_idx, ampcopy in enumerate(sorted(amp_wins_with_masked.keys())):
				plt.subplot(copy_numbers[amplicon], 1, ampcopy_idx + 1)
				
				plt.plot([x for x in xrange(len(amp_wins_with_masked[ampcopy])) if max_yval > amp_wins_with_masked[ampcopy][x] >=0], 
						 [x for x in amp_wins_with_masked[ampcopy] if max_yval > x >= 0], 
						 'b.', markersize=1)
				plt.plot([x for x in xrange(len(amp_wins_with_masked[ampcopy])) if amp_wins_with_masked[ampcopy][x] >= max_yval], 
				[max_yval for x in amp_wins_with_masked[ampcopy] if x >= max_yval], 
				'cD', markersize=1)
			
		
				if args.segmentation:
				
					breakpoint = maxps[ampcopy][1]
					plt.plot([0, breakpoint], [maxps[ampcopy][2], maxps[ampcopy][2]], 'r', lw=.4)
					plt.plot([breakpoint, len(amp_wins_with_masked[ampcopy])], [maxps[ampcopy][3], maxps[ampcopy][3]], 'r', lw=.4)
					plt.plot([breakpoint, breakpoint], [0, max_yval], 'r', lw=.9)

					plt.text(.01 * len(amp_wins_with_masked[ampcopy]), .96 * max_yval, 
							 str(int(maxps[ampcopy][4])), size=6, 
							 horizontalalignment='left', 
							 verticalalignment='top')
					plt.text(.99 * len(amp_wins_with_masked[ampcopy]), .96 * max_yval, 
							 str(int(maxps[ampcopy][5])), size=6, 
							 horizontalalignment='right', 
							 verticalalignment='top')

				
				plt.axis([0, len(amp_wins_with_masked[ampcopy]), 0, max_yval])
				
				if ampcopy_idx != copy_numbers[amplicon] - 1:
					plt.xticks([])
			
 				if args.segmentation:
 					plt.title('%s   %.2e' %(ampcopy, maxps[ampcopy][0]), size=7, y=.96)
 				else:
					plt.title('%s' %(ampcopy), size=7, y=.96)
			
			if args.segmentation:
				plt.savefig('%s_%s_segmented.png' %(args.outfile, amplicon), dpi=500)
			else:
				plt.savefig('%s_%s.png' %(args.outfile, amplicon), dpi=500)
			
