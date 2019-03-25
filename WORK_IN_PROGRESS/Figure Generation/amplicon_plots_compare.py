#!/usr/bin/python

'''
Generate chrY amplicon CNV comparison plot
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
others = ['TSPY', 'DYZ19', 'RBMY', 'DAZ']
controls = ['Ctrl_reg_0', 'Ctrl_reg_1', 'Ctrl_reg_2', 'Ctrl_reg_3', 'Ctrl_reg_4', 'Ctrl_reg_5']


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

#For new method of calling "wrong" controls. These are the means and stds of the four large control
#regions in the 1000 Genomes data.
all_file = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons.pickle')
ctrl_means = {'Ctrl_reg_0': np.mean(all_file['Ctrl_reg_0']), 
			  'Ctrl_reg_1': np.mean(all_file['Ctrl_reg_1']), 
			  'Ctrl_reg_2': np.mean(all_file['Ctrl_reg_2']), 
			  'Ctrl_reg_5': np.mean(all_file['Ctrl_reg_5'])}
ctrl_stds = {'Ctrl_reg_0': np.std(all_file['Ctrl_reg_0']), 
			 'Ctrl_reg_1': np.std(all_file['Ctrl_reg_1']), 
			 'Ctrl_reg_2': np.std(all_file['Ctrl_reg_2']), 
			 'Ctrl_reg_5': np.std(all_file['Ctrl_reg_5'])}
			  
			  
def hexcolor_adjust(hexcolor, adjustment):
	#Adjustment <1 darkens, >1 lightens
	ints = [int(hexcolor[1:3],16), int(hexcolor[3:5],16), int(hexcolor[5:7],16)]
	adjints = [min(255,int(x*adjustment)) for x in ints]
	adjhex = [hex(x)[2:] if len(hex(x)) == 4 else '0' + hex(x)[2:] for x in adjints]
	return '#%s%s%s' %(adjhex[0], adjhex[1], adjhex[2])


def plot_man_compare(amp_vals1, amp_vals2, man_id, plot_all=False):
	plt.figure(figsize=[5.77, 5.64])

	if not plot_all:
		amps = inv_reps + ['gap'] + pals + ['gap'] + azfc + ['gap']
	else:
		amps = inv_reps + ['gap'] + pals + ['gap'] + azfc + ['gap'] + others + ['gap'] + controls + ['gap']

	mid_offset = 1
	if not plot_all:
		width = .25
	else:
		width = .15
	ax = plt.subplot(polar=True)
	ax.axis('Off')
	gray_color = '#EEEEEE'
	
	#Make background gray circle
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
				size=5.9 if not plot_all else 4, 
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
		bar_height1 = amp_vals1[amp_idx]
		bar_color = '#4C7A8A'
		if amp in inv_reps + pals + azfc:
			if int(bar_height1 * copy_numbers[amp] + 0.5) != copy_numbers[amp]:
				bar_color = '#F67171' if bar_height1 < 1 else '#88E26B'
		if amp in ['Ctrl_reg_0', 'Ctrl_reg_1', 'Ctrl_reg_2', 'Ctrl_reg_5'] and abs(bar_height1 - ctrl_means[amp]) > 2 * ctrl_stds[amp]:
			bar_color = '#F67171'
		bar_color = hexcolor_adjust(bar_color, 1.4)
		ax.fill([bar_center - np.pi/2,
					bar_center - np.arctan((width/2) / (bar_height1 + mid_offset)),
					bar_center,
					bar_center],
				[width/2,
					np.sqrt((width/2)**2 + (bar_height1 + mid_offset)**2), 
					bar_height1 + mid_offset, 
					0],
				fc = bar_color,
				lw = 0,
				zorder = 3)
		
		bar_height2 = amp_vals2[amp_idx]
		bar_color = '#4C7A8A'
		if amp in inv_reps + pals + azfc:
			if int(bar_height2 * copy_numbers[amp] + 0.5) != copy_numbers[amp]:
				bar_color = '#F67171' if bar_height2 < 1 else '#88E26B'
		if amp in ['Ctrl_reg_0', 'Ctrl_reg_1', 'Ctrl_reg_2', 'Ctrl_reg_5'] and abs(bar_height2 - ctrl_means[amp]) > 2 * ctrl_stds[amp]:
			bar_color = '#F67171'
		bar_color = hexcolor_adjust(bar_color, 0.6)
		ax.fill([bar_center,
					bar_center,
					bar_center + np.arctan((width/2) / (bar_height2 + mid_offset)),
					bar_center + np.pi/2],
				[0,
					bar_height2 + mid_offset, 
					np.sqrt((width/2)**2 + (bar_height2 + mid_offset)**2), 
					width/2],
				fc = bar_color,
				lw = 0,
				zorder = 3)
				
				
	#Make lines representing copy number
	inner_offset = 0.08
	for amp_idx, amp in enumerate(amps):
		if amp == 'gap' or amp in others + controls:
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
		if amp_vals1[amp_idx] < 2 and amp_vals2[amp_idx] < 2:
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
	
	if os.path.isfile(sys.argv[3]):
		print 'OUTFILE ALREADY EXISTS'
		#quit()
		
	plot_all = False
	if not plot_all:
		amps = inv_reps + ['gap'] + pals + ['gap'] + azfc + ['gap']
	else:
		amps = inv_reps + ['gap'] + pals + ['gap'] + azfc + ['gap'] + others + ['gap'] + controls + ['gap']

	plt.clf()
	plt.figure()
	man1 = sys.argv[1]
	man2 = sys.argv[2]
	vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_amplicons.pickle')
	if man1 in vectors.index:
		man1_vec = vectors.loc[man1][amps]
		man1_name = man1
	else:
		man1_vec = vecfile_to_vector(man1).loc[0][amps]
		#man1_name = '.'.join(man1.split('/')[-1].split('.')[:-1])
		man1_name = man1.split('/')[-1].split('_')[0]
		
	if man2 in vectors.index:
		man2_vec = vectors.loc[man2][amps]
		man2_name = man2
	else:
		man2_vec = vecfile_to_vector(man2).loc[0][amps]
		#man2_name = '.'.join(man2.split('/')[-1].split('.')[:-1])
		man2_name = man2.split('/')[-1].split('_')[0]
		
	#plot_man_compare(man1_vec, man2_vec, '%s & %s' %(man1_name, man2_name), plot_all=plot_all)
	#plot_man_compare(man1_vec, man2_vec, sys.argv[1].split('/')[0], plot_all=True)
	#plot_man_compare(man1_vec, man2_vec, sys.argv[1].split('/')[-1][:7], plot_all=True)
	plot_man_compare(man1_vec, man2_vec, man1_name, plot_all=False)
	
	plt.savefig(sys.argv[3], dpi=500, bbox_inches='tight')
		




