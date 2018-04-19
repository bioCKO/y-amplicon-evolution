#!/usr/bin/python

'''
Make figures of most significant CBS hits
'''

import pickle
import glob
import sys
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from scipy.stats import mannwhitneyu
#from make_chrY_vector import get_reg_depth
#from amplicon_CBS import bed_to_depthvec, get_reg_wins, get_amplicon_wins


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



def multibed_to_depthvec(*infiles):
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
	return chrY_depth


def get_loc_corrs(winstart, winend, mask_file):
	#Get conversion locations from skipping masked to not skipping masked
	loc_corrs = {}
	noskipwins = get_reg_wins([0] * 27000000, winstart, winend, mask_repeats=mask_file, skip_masked=False)
	skipdex = 0
	for windex, win in enumerate(noskipwins):
		if win != -1:
			loc_corrs[skipdex] = windex
			skipdex += 1
	return loc_corrs
	
def get_amp_corrs(rep_dict, mask_file, *ampregs):
	#Get conversion locations from skipping masked to not skipping masked for multi_region amplicons
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
	
	out_dir = '/lab/Page_lab-users/lsteitz/Amplicon_Edge_Figures'
	
	goodlocsfile = open('/lab/solexa_page/lsteitz/Y_repetition/good_chrY_amplicon_bases_w100_c11.pickle', 'r')
	goodlocs = pickle.load(goodlocsfile)
	goodlocsfile.close()
	
	cbsfile = open('/lab/solexa_page/lsteitz/1000_Ys/partial_cnvs/CBS_analysis_ext_all.txt', 'r')
	cbs = []
	for line in cbsfile:
		cbs.append(line.rstrip().split('\t'))
	cbsfile.close()
	
	zzz = 'HG02684    Green   g3      0.0000e+00      1711    2.878   7.120   2       6       g2      0.0000e+00      1808    3.049   7.799   2       6       g1      0.0000e+00      1673    2.943   7.404   2       6'
	zzz = zzz.split()
	cbs = [zzz]
	
	regs = {x: [] for x in list(set([x[1] for x in cbs]))}
	
	for reg in regs:
		regcbs = [x for x in cbs if x[1] == reg]
		if reg.startswith('Ctrl') or reg.startswith('Norm'):
			regps = [[x[0], float(x[2]), x] for x in regcbs]
		else:
			regps = [[x[0], max([float(x[y]) for y in xrange(3, len(x), 7)]), x] for x in regcbs]
		regs[reg] = regps


	lowmen = {}

	for reg in regs:
		lowmen[reg] = sorted(regs[reg], key=lambda x: x[1])[:64]
	
	
	ctrl_regs = {'Ctrl_reg_0': [4000000, 5000000], 
				 'Ctrl_reg_1': [7000000, 7500000], 
				 'Ctrl_reg_2': [8000000, 8500000], 
				 'Ctrl_reg_3': [16780000, 16800000], 
				 'Ctrl_reg_4': [19550000, 19560000], 
				 'Ctrl_reg_5': [20500000, 20550000]}
	
	
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
	
	
	if sys.argv[1] == '-s':
		plot_splits = True
		man = sys.argv[2]
		regs = sys.argv[3:]
	else:
		plot_splits = False
		man = sys.argv[1]
		regs = sys.argv[2:]
		
	for reg in regs:
		
		plt.clf()
		
		#THE NORM AND CTRL REG CODE IS FROM CBS_top_hits_figs_v2.py; IT WILL NOT WORK IN THIS SCRIPT
		#IF THE CONTINUE STATEMENTS ARE DELETED
		if reg == 'Norm reg':
			
			###
			continue
			###
			
			loc_corrs = get_loc_corrs(14500000, 15500000, goodlocs)
			
			plt.figure(figsize=[30, 24])
			for fig_idx, man in enumerate(lowmen[reg]):
				plt.subplot(16, 4, fig_idx + 1)

				infiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/*/%s/*/*_chrY_cov_rl50plus_pe.gccorrected.bedgraph' %(man[0]))		
				chrY_depth = multibed_to_depthvec(*infiles)
			

				normwins = get_reg_wins(chrY_depth, 14500000, 15500000, mask_repeats=goodlocs, skip_masked=False)#, winsize=10000)
				
				max_yval = int(max(float(man[2][4]), float(man[2][5])) + 3 * np.std(normwins)) + 1

				plt.plot([x for x in xrange(len(normwins)) if max_yval > normwins[x] >=0], [x for x in normwins if max_yval > x >= 0], 'b.', markersize=1)
				plt.plot([x for x in xrange(len(normwins)) if normwins[x] >= max_yval], [max_yval for x in normwins if x >= max_yval], 'cD', markersize=1)
				
				breakpoint = loc_corrs[int(man[2][3])]
				
				plt.plot([0, breakpoint], [float(man[2][4]), float(man[2][4])], 'r', lw=.4)
				plt.plot([breakpoint, len(normwins)], [float(man[2][5]), float(man[2][5])], 'r', lw=.4)
				plt.plot([breakpoint, breakpoint], [0, max_yval], 'r', lw=.9)
				
				plt.axis([0, len(normwins), 0, max_yval])
				if fig_idx < 61:
					plt.xticks([])
				plt.title('%s %.2e' %(man[0], man[1]), size=7, y=.96)

			plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/Norm_wins.png', dpi=500)
		
		
		elif reg.startswith('Ctrl'):
		
			###
			continue
			###
			
			
			cstart = ctrl_regs[reg][0]
			cend = ctrl_regs[reg][1]
			loc_corrs = get_loc_corrs(cstart, cend, goodlocs)
			
			plt.figure(figsize=[30, 24])
			for fig_idx, man in enumerate(lowmen[reg]):
				plt.subplot(16, 4, fig_idx + 1)

				infiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/*/%s/*/*_chrY_cov_rl50plus_pe.gccorrected.bedgraph' %(man[0]))		
				chrY_depth = multibed_to_depthvec(*infiles)
		
				ctrlwins = get_reg_wins(chrY_depth, cstart, cend, mask_repeats=goodlocs, skip_masked=False)
				max_yval = int(max(float(man[2][4]), float(man[2][5])) + 3 * np.std(ctrlwins)) + 1
		
				plt.plot([x for x in xrange(len(ctrlwins)) if max_yval > ctrlwins[x] >=0], [x for x in ctrlwins if max_yval > x >= 0], 'b.', markersize=1)
				plt.plot([x for x in xrange(len(ctrlwins)) if ctrlwins[x] >= max_yval], [max_yval for x in ctrlwins if x >= max_yval], 'cD', markersize=1)
		
				breakpoint = loc_corrs[int(man[2][3])]
				
				plt.plot([0, breakpoint], [float(man[2][4]), float(man[2][4])], 'r', lw=.4)
				plt.plot([breakpoint, len(ctrlwins)], [float(man[2][5]), float(man[2][5])], 'r', lw=.4)
				plt.plot([breakpoint, breakpoint], [0, max_yval], 'r', lw=.9)
				
				plt.axis([0, len(ctrlwins), 0, max_yval])
				if fig_idx < 61:
					plt.xticks([])
				plt.title('%s %.2e' %(man[0], man[1]), size=7, y=.96)

			plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s.png' %(reg), dpi=500)
		
		
		else:
			loc_corrs = get_amp_corrs(y_reps, goodlocs, *amplicon_regions[reg])
			
			plot_count = 1
			
			plt.figure(figsize=[7.5, 1.5 * copy_numbers[reg]])
			infiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/*/%s/*/*_chrY_cov_rl50plus_pe.gccorrected.bedgraph' %(man))
			#infiles = glob.glob('/lab/solexa_page/lsteitz/A00/%s.chrY_to_hg38.bedgraph' %(man))
			chrY_depth = multibed_to_depthvec(*infiles)
			ampwins = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[reg], mask_repeats=goodlocs, skip_masked=False)
			
			max_yval = int(max([np.mean(ampwins[x]) for x in ampwins]) + 3 * max([np.std(ampwins[x]) for x in ampwins]) + 1)
			
			for ampcopy_idx, ampcopy in enumerate(sorted(ampwins.keys())):
				plt.subplot(copy_numbers[reg], 1, ampcopy_idx + 1)
				
				plt.plot([x for x in xrange(len(ampwins[ampcopy])) if max_yval > ampwins[ampcopy][x] >=0], 
						 [x for x in ampwins[ampcopy] if max_yval > x >= 0], 
						 'b.', markersize=1)
				plt.plot([x for x in xrange(len(ampwins[ampcopy])) if ampwins[ampcopy][x] >= max_yval], 
				[max_yval for x in ampwins[ampcopy] if x >= max_yval], 
				'cD', markersize=1)
				
				
				if plot_splits:
				
					man_data_raw = next(x for x in lowmen[reg] if x[0] == man)
					man_data = {man_data_raw[2][x]: map(float, man_data_raw[2][x+1:x+7]) for x in xrange(2,len(man_data_raw[2]),7)}

					breakpoint = loc_corrs[ampcopy][man_data[ampcopy][1]]
					
					plt.plot([0, breakpoint], [man_data[ampcopy][2], man_data[ampcopy][2]], 'r', lw=.4)
					plt.plot([breakpoint, len(ampwins[ampcopy])], [man_data[ampcopy][3], man_data[ampcopy][3]], 'r', lw=.4)
					plt.plot([breakpoint, breakpoint], [0, max_yval], 'r', lw=.9)
				
					plt.text(.01 * len(ampwins[ampcopy]), .96 * max_yval, 
							 str(int(man_data[ampcopy][4])), size=6, 
							 horizontalalignment='left', 
							 verticalalignment='top')
					plt.text(.99 * len(ampwins[ampcopy]), .96 * max_yval, 
							 str(int(man_data[ampcopy][5])), size=6, 
							 horizontalalignment='right', 
							 verticalalignment='top')

			
				plt.axis([0, len(ampwins[ampcopy]), 0, max_yval])
				
				if ampcopy_idx != copy_numbers[reg] - 1:
					plt.xticks([])
				
				if ampcopy_idx == 0:	
					if plot_splits:
						plt.title('%s\n%s   %.2e' %(man, ampcopy, man_data[ampcopy][0]), size=7, y=.96)
					else:
						plt.title('%s\n%s' %(man, ampcopy), size=7, y=.96)
				
				else:
					if plot_splits:
						plt.title('%s   %.2e' %(ampcopy, man_data[ampcopy][0]), size=7, y=.96)
					else:
						plt.title('%s' %(ampcopy), size=7, y=.96)
			
			if plot_splits:
				#plt.savefig('%s/%s_%s_split.svg' %(out_dir, man, reg))
				plt.savefig('%s/%s_%s_split.png' %(out_dir, man, reg), dpi=500)
			else:
				#plt.savefig('%s/%s_%s.svg' %(out_dir, man, reg), dpi=500)
				plt.savefig('%s/%s_%s.png' %(out_dir, man, reg), dpi=500)
				
				
				
				
				
			
			'''
			
			for man_idx, man in enumerate(lowmen[reg]):
				
				man_data = {man[2][x]: map(float, man[2][x+1:x+7]) for x in xrange(2,len(man[2]),7)}
				
				
				men_per_fig = (16 / (len(man_data) + 1)) * 4
				if man_idx%men_per_fig == 0 and man_idx > 0:
					plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s_%i.png' %(reg, plot_count), dpi=500)
					plot_count += 1
					plt.clf()
					plt.figure(figsize=[30, 24])
									
				
				infiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/*/%s/*/*_chrY_cov_rl50plus_pe.gccorrected.bedgraph' %(man[0]))		
				chrY_depth = multibed_to_depthvec(*infiles)
				
				ampwins = get_amplicon_wins(chrY_depth, y_reps, *amplicon_regions[reg], mask_repeats=goodlocs, skip_masked=False)
				max_yval = int(max([max(man_data[x][2:4]) for x in man_data]) + 3 * max([np.std(ampwins[x]) for x in ampwins]) + 1)
								
				for ampcopy_idx, ampcopy in enumerate(sorted(man_data.keys())):
				
					subplot_idx = 1 + man_idx%men_per_fig + (4 * ampcopy_idx) + (len(man_data) * 4 *((man_idx%men_per_fig)/4))
					plt.subplot(16, 4, subplot_idx)
					
					plt.plot([x for x in xrange(len(ampwins[ampcopy])) if max_yval > ampwins[ampcopy][x] >=0], 
							 [x for x in ampwins[ampcopy] if max_yval > x >= 0], 
							 'b.', markersize=1)
					plt.plot([x for x in xrange(len(ampwins[ampcopy])) if ampwins[ampcopy][x] >= max_yval], 
					[max_yval for x in ampwins[ampcopy] if x >= max_yval], 
					'cD', markersize=1)
		
					breakpoint = loc_corrs[ampcopy][man_data[ampcopy][1]]
					
					plt.plot([0, breakpoint], [man_data[ampcopy][2], man_data[ampcopy][2]], 'r', lw=.4)
					plt.plot([breakpoint, len(ampwins[ampcopy])], [man_data[ampcopy][3], man_data[ampcopy][3]], 'r', lw=.4)
					plt.plot([breakpoint, breakpoint], [0, max_yval], 'r', lw=.9)
					
					plt.text(.01 * len(ampwins[ampcopy]), .96 * max_yval, 
							 str(int(man_data[ampcopy][4])), size=6, 
							 horizontalalignment='left', 
							 verticalalignment='top')
					plt.text(.99 * len(ampwins[ampcopy]), .96 * max_yval, 
							 str(int(man_data[ampcopy][5])), size=6, 
							 horizontalalignment='right', 
							 verticalalignment='top')
				
					plt.axis([0, len(ampwins[ampcopy]), 0, max_yval])
					
					if ampcopy_idx != len(man_data) - 1:
						plt.xticks([])
					
					if ampcopy_idx == 0:	
						#plt.title(man[0], size=12)
						plt.title('%s\n%s   %.2e' %(man[0], ampcopy, man_data[ampcopy][0]), size=7, y=.96)
					else:
						plt.title('%s   %.2e' %(ampcopy, man_data[ampcopy][0]), size=7, y=.96)					
						
			plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Partial_SVs/%s_%i.png' %(reg, plot_count), dpi=500)
			
			'''
	
	
	
	
	






















