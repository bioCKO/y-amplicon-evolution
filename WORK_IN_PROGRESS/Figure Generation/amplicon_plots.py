#!/usr/bin/python

'''
Generate chrY amplicon CNV plot
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import numpy as np
import sys
import pandas as pd
import os
import math


inv_reps = ['IR1', 'IR2', 'IR3', 'IR5']
pals = ['P4', 'P5', 'P6', 'P7', 'P8']
azfc = ['Blue', 'Teal', 'Green', 'Red', 'Gray', 'Yellow']
amps = inv_reps + ['gap'] + pals + ['gap'] + azfc + ['gap']

cn_2 = ['IR3', 'IR1', 'P8', 'P7', 'P6', 'P5', 'P4', 'IR2', 'Teal', 'Gray', 'Yellow']
cn_3 = ['Green']
cn_4 = ['IR5', 'Blue', 'Red']

copy_numbers = {}
for amp in cn_2:
	copy_numbers[amp] = 2
for amp in cn_3:
	copy_numbers[amp] = 3
for amp in cn_4:
	copy_numbers[amp] = 4

max_cov = 3.5

def plot_man(amp_vals, man_id):
	plt.figure(figsize=[5.77, 5.64])
	mid_offset = 1
	width = .25
	ax = plt.subplot(polar=True)
	ax.axis('Off')
	
	
	
	#Make background gray circle
	gray_color = '#EEEEEE'
	ax.bar(0 - np.pi/len(amps), 
		   2, 
		   bottom = mid_offset, 
		   width = 2*np.pi, 
		   fc=gray_color, 
		   lw=0, 
		   zorder=1)
	
	#Make central background white circle
	ax.bar(0 - np.pi/len(amps), 
		   mid_offset, 
		   width = 2*np.pi, 
		   fc='white', 
		   lw=0, 
		   zorder=4)
	
	#Make circles at 0x and 1x coverage and inside text location
	circle_rs = np.arange(0,2*np.pi,.001) - (np.pi/len(amps)) #The last term is so that in non-polar mode, the lines start where the wedges start
	plt.plot(circle_rs, [mid_offset]*len(circle_rs), color='black', zorder=5, lw=1)
	plt.plot(circle_rs, [mid_offset + 1]*len(circle_rs), color='gray', ls='dashed', zorder=2, lw=1, dashes=(6,4.9))
	
	#Make labels of amplicon names within inner circle
	text_locs = np.arange(0, 2*np.pi, 2*np.pi/len(amps))
	text_mid_offset = 0.04
	for amp_idx, amp in enumerate(amps):
		if amp == 'gap':
			continue
		ax.text(text_locs[amp_idx], 
				mid_offset - text_mid_offset, 
				amp, 
				size=5.9, 
				rotation= 360. / len(amps) * amp_idx, 
				va='center', 
				ha='right', 
				rotation_mode='anchor',
				zorder=5)		
	
	
	
	#Locations of bar centers around the plot
	bar_locs = np.arange(0, 2*np.pi, 2*np.pi/len(amps))
	
	
	#Make bars
	for amp_idx, amp in enumerate(amps):
		if amp == 'gap':
			continue
		bar_center = bar_locs[amp_idx]
		bar_height = amp_vals[amp_idx]
		if int(amp_vals[amp_idx] * copy_numbers[amp] + 0.5) != copy_numbers[amp]:
			bar_color = '#F67171' if amp_vals[amp_idx] < 1 else '#88E26B'
		else:
			bar_color = '#4C7A8A'
		ax.fill([bar_center - np.pi/2,
					bar_center - np.arctan((width/2) / (bar_height + mid_offset)),
					bar_center + np.arctan((width/2) / (bar_height + mid_offset)),
					bar_center + np.pi/2],
				[width/2,
					np.sqrt((width/2)**2 + (bar_height + mid_offset)**2), 
					np.sqrt((width/2)**2 + (bar_height + mid_offset)**2), 
					width/2],
				fc = bar_color,
				lw = 0,
				zorder = 3)
				
				
	#Make lines representing copy number
	inner_offset = 0.08
	for amp_idx, amp in enumerate(amps):
		if amp == 'gap':
			continue
		for cn in range(1,copy_numbers[amp] * 2 + 1):
			line_center = bar_locs[amp_idx]
			line_height = mid_offset + (1./copy_numbers[amp]) * cn - (1./copy_numbers[amp])/2
			ax.plot([line_center - np.arctan(((width - inner_offset) / 2) / line_height),
					line_center + np.arctan(((width - inner_offset) / 2) / \
					line_height)],
					[np.sqrt(((width - inner_offset) / 2)**2 + line_height**2)] * 2, 
					color = gray_color, lw=1.5, zorder=4)
		#Lines outside gray circle, for very high copy numbers:	
		if amp_vals[amp_idx] < 2:
			continue
		for cn in range(copy_numbers[amp] * 2 + 1, copy_numbers[amp] * 4):
			line_center = bar_locs[amp_idx]
			line_height = mid_offset + (1./copy_numbers[amp]) * cn - (1./copy_numbers[amp])/2
			ax.plot([line_center - np.arctan(((width - inner_offset) / 2) / line_height),
					line_center + np.arctan(((width - inner_offset) / 2) / \
					line_height)],
					[np.sqrt(((width - inner_offset) / 2)**2 + line_height**2)] * 2, 
					color = '#FFFFFF', lw=1.5, zorder=4)
	
	#Write individual ID
	ax.text(3 * np.pi / 2, 
			4.3, 
			man_id,
			size=10,
			ha='center',
			zorder=5)	
	
	ax.set_ylim(0, max_cov + mid_offset)

def vecfile_to_vector(vecfile):
	with open(vecfile, 'r') as vecin:
		veclines = vecin.readlines()
	vals = [map(float, veclines[2].rstrip().split())]
	names = veclines[1][1:].rstrip().split()
	return pd.DataFrame(vals, columns=names)

	

if __name__ == '__main__':
	
	#man = sys.argv[1]
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls.pickle')
	vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons.pickle')
	
	
	men = sys.argv[1:-1]
	for man in men:
		if os.path.isfile('%s%s.png' %(sys.argv[-1], man)):
			print 'OUTFILE ALREADY EXISTS'
		if man in vectors.index:
			man_vec = vectors.loc[man][amps]
			plot_man(man_vec, man)
			plt.savefig('%s%s.svg' %(sys.argv[-1], man), dpi=500, bbox_inches='tight')
		else:
			man_vec = vecfile_to_vector(man).loc[0][amps]
			#man_name = '.'.join(sys.argv[-1].split('/')[-1].split('.')[:-1])
			man_name = man.split('/')[-1].split('_')[0]
			plot_man(man_vec, man_name)
			plt.savefig('%s%s.svg' %(sys.argv[-1], man_name), dpi=500, bbox_inches='tight')
			
		#out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Noncontrol Men Circle Plots/'
		#plt.savefig('%s%s.png' %(out_dir, sys.argv[1]), dpi=500, bbox_inches='tight')
		#plt.savefig('%s%s.png' %(sys.argv[-1], man), dpi=500, bbox_inches='tight')
		#plt.savefig('%s.png' %(sys.argv[-1]), dpi=500, bbox_inches='tight')
		#plt.savefig('%s.svg' %(sys.argv[-1]), bbox_inches='tight')
		




