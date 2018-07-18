#!/usr/bin/python

'''
Several analyses on copy number calls:
Check that lower depth isn't correlated with increased CNV calls
Get percantage of deletions, duplications, and complex variants in the whole dataset and in each haplogroup
'''

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import numpy as np

cn_2 = ['IR3', 'IR1', 'P8', 'P7', 'P6', 'P5', 'P4', 'IR2', 'Teal', 'Gray', 'Yellow']
cn_3 = ['Green']
cn_4 = ['IR5', 'Blue', 'Red']
complex_and_control = ['TSPY', 'DYZ19', 'RBMY', 'DAZ', #'P7',
					   'Ctrl_reg_0', 'Ctrl_reg_1', 
					   'Ctrl_reg_2', 'Ctrl_reg_3', 'Ctrl_reg_4', 
					   'Ctrl_reg_5']#, 'IR3', 'IR1', 'P8', 'P6', 'P5', 'P4', 'IR2', 'IR5']

copy_numbers = {}
for amp in cn_2:
	copy_numbers[amp] = 2
for amp in cn_3:
	copy_numbers[amp] = 3
for amp in cn_4:
	copy_numbers[amp] = 4


vectors = pd.read_pickle(sys.argv[1])
del vectors['P7']


#Check that lower depth doesn't cause more called SVs
#Not controlling for differences of frequency within haplogroups
#Looking at binary: if at least 1 SV called

from statsmodels.stats.proportion import proportion_confint

sv_cutoffs = [2, 3, 4, 5, 10000] #This list NEEDS TO BE SORTED FROM LOWEST TO HIGHEST!
sv_cutoffs.sort() #Just to make sure
sv_counter = {x:[0,0] for x in sv_cutoffs}
#Key values are coverage cutoffs (the last one should be higher than any possible coverage); 
#value lists are [total, number with at least one SV]

for ind in vectors.index:
	depth = vectors.loc[ind, 'Single-copy_depth']
	for cutoff in sv_cutoffs:
		if depth < cutoff:
			sv_counter[cutoff][0] += 1
			break
	for amp in vectors.columns[3:]:
		if vectors.loc[ind, amp] != copy_numbers[amp]:
			for cutoff in sv_cutoffs:
				if depth < cutoff:
					sv_counter[cutoff][1] += 1
					break
			break

sv_calls = sum([x[1] for x in sv_counter.values()])
print 'Number of individuals with at least one SV: %i out of %i (%i%%)' %(sv_calls, len(vectors), float(sv_calls)/len(vectors)*100)

last_cutoff = 0
for cutoff in sv_cutoffs:
	sv_percent = int(sv_counter[cutoff][1]/float(sv_counter[cutoff][0]) * 100)
	print 'Percent of individuals with %i < coverage < %i with SVs: %i%% (%i out of %i)' %(last_cutoff, cutoff, sv_percent, sv_counter[cutoff][1], sv_counter[cutoff][0])
	last_cutoff = cutoff
print

def hexcolor_adjust(hexcolor, adjustment):
	#Adjustment <1 darkens, >1 lightens
	ints = [int(hexcolor[1:3],16), int(hexcolor[3:5],16), int(hexcolor[5:7],16)]
	adjints = [min(255,int(x*adjustment)) for x in ints]
	adjhex = [hex(x)[2:] if len(hex(x)) == 4 else '0' + hex(x)[2:] for x in adjints]
	return '#%s%s%s' %(adjhex[0], adjhex[1], adjhex[2])

#adjs = [.4, .8, 1.2, 1.6, 2]
adjs = [.4, .7, 1, 1.3, 1.6]
graycode = '#808080'

plt.clf()
fig = plt.figure(figsize=[4.8, 6.4])
ax = fig.add_subplot(111)
last_cutoff = 0
for cutidx, cutoff in enumerate(sv_cutoffs):
	y_error = proportion_confint(sv_counter[cutoff][1], sv_counter[cutoff][0])
	sv_pct = sv_counter[cutoff][1]/float(sv_counter[cutoff][0])
	ax.bar(cutidx+0.5, sv_pct, color=hexcolor_adjust(graycode, adjs[cutidx]),
			yerr=[[sv_pct-y_error[0]], [y_error[1]-sv_pct]], ecolor='black')
ax.xaxis.set_tick_params(width=0)
plt.xticks(np.arange(.5,6,1), ['0x-2x', '2x-3x', '3x-4x', '4x-5x', '5x+'], size=12)
plt.yticks(size=12)
plt.ylabel('% individuals with CNV', size=16)
plt.xlabel('Coverage', size=16)
plt.axis([-0.2,5.2,0,.3])
#plt.axis([-0.2,5.2,0,.8])

plt.tick_params(axis='x', direction='out', top='off')
plt.tick_params(axis='y', direction='out', right='off')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
	
plt.savefig('./depth_control.svg', bbox_inches='tight')
plt.savefig('./depth_control.png', dpi=300, bbox_inches='tight')


#Check within large haplogroups to control for differing frequencies of SV between haplogroups:

haps = ['R1a', 'R1b', 'O3', 'E1b']
for hap in haps:
	sv_counter = {x:[0,0] for x in sv_cutoffs}
	#Key values are coverage cutoffs (the last one should be higher than any possible coverage); 
	#value lists are [total, number with at least one SV]

	for ind in vectors.index:
		if vectors.loc[ind, 'Sub-Haplogroup'] != hap:
			continue
		depth = vectors.loc[ind, 'Single-copy_depth']
		for cutoff in sv_cutoffs:
			if depth < cutoff:
				sv_counter[cutoff][0] += 1
				break
		for amp in vectors.columns[3:]:
			if vectors.loc[ind, amp] != copy_numbers[amp]:
				for cutoff in sv_cutoffs:
					if depth < cutoff:
						sv_counter[cutoff][1] += 1
						break
				break

	sv_calls = sum([x[1] for x in sv_counter.values()])
	print 'Number of individuals in %s with at least one SV: %i out of %i (%i%%)' \
			%(hap, 
				sv_calls, 
				len(vectors.loc[vectors['Sub-Haplogroup'] == hap]), 
				float(sv_calls)/len(vectors.loc[vectors['Sub-Haplogroup'] == hap])*100
				)

	last_cutoff = 0
	for cutoff in sv_cutoffs:
		if sv_counter[cutoff][0] == 0:
			print 'Percent of individuals in %s with %i < coverage < %i with SVs: N/A (0 out of 0)' \
				%(hap, last_cutoff, cutoff)
			continue
		sv_percent = int(sv_counter[cutoff][1]/float(sv_counter[cutoff][0]) * 100)
		print 'Percent of individuals in %s with %i < coverage < %i with SVs: %i%% (%i out of %i)' \
			%(hap, last_cutoff, cutoff, sv_percent, sv_counter[cutoff][1], sv_counter[cutoff][0])
		last_cutoff = cutoff
	print

	hapcolors = {'A': {'color': '#F79981', 'haplogroups': ['A0', 'A1a']},
		'C': {'color': '#35BEFC', 'haplogroups': ['C5', 'C3', 'C1']},
		'B': {'color': '#6C670A', 'haplogroups': ['B2', 'B']},
		'E': {'color': '#FF19B0', 'haplogroups': ['E1b', 'E2', 'E1a']},
		'D': {'color': '#073018', 'haplogroups': ['D2']},
		'G': {'color': '#F7FB92', 'haplogroups': ['G1', 'G2']},
		'F': {'color': '#6EFFB6', 'haplogroups': ['F']},
		'I': {'color': '#23233B', 'haplogroups': ['I2', 'I1']},
		'H': {'color': '#0276A4', 'haplogroups': ['H', 'H1', 'H0', 'H2']},
		'J': {'color': '#900A2A', 'haplogroups': ['J1', 'J2']},
		'L': {'color': '#FFDCDE', 'haplogroups': ['L1']},
		'O': {'color': '#6254C5', 'haplogroups': ['O3', 'O2', 'O1']},
		'N': {'color': '#C31929', 'haplogroups': ['NO', 'N', 'N1']},
		'Q': {'color': '#4AF3F2', 'haplogroups': ['Q1b', 'Q1a']},
		'R': {'color': '#A56410', 'haplogroups': ['R1a', 'R1b', 'R2']},
		'T': {'color': '#BDFD40', 'haplogroups': ['T']},
		'NONE': {'color': '#024444', 'haplogroups': ['NONE']}}
	
	adjs = [.4, .8, 1.2, 1.6, 2]
	colorcode = hapcolors[hap[0]]['color']
	
	plt.clf()
	fig = plt.figure(figsize=[4.8, 6.4])
	ax = fig.add_subplot(111)
	last_cutoff = 0
	for cutidx, cutoff in enumerate(sv_cutoffs):
		y_error = proportion_confint(sv_counter[cutoff][1], sv_counter[cutoff][0])
		sv_pct = sv_counter[cutoff][1]/float(sv_counter[cutoff][0])
		ax.bar(cutidx+0.5, sv_pct, color=hexcolor_adjust(colorcode, adjs[cutidx]),
				yerr=[[min(sv_pct,sv_pct-y_error[0])], [y_error[1]-sv_pct]], ecolor='black')
	ax.xaxis.set_tick_params(width=0)
	plt.xticks(np.arange(.5,6,1), ['0x-2x', '2x-3x', '3x-4x', '4x-5x', '5x+'], size=12)
	plt.ylabel('%% %s individuals with CNV' %(hap), size=16)
	plt.xlabel('Coverage', size=16)
	plt.yticks(size=12)
	plt.axis([-0.2,5.2,0,0.3])
	#plt.axis([-0.2,5.2,0,0.8])

	plt.tick_params(axis='x', direction='out', top='off')
	plt.tick_params(axis='y', direction='out', right='off')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	plt.savefig('./depth_control_%s.svg' %(hap), bbox_inches='tight')
	plt.savefig('./depth_control_%s.png' %(hap), dpi=300, bbox_inches='tight')


#See if deletions or duplications are more common

normal_cn = [copy_numbers[x] for x in vectors.columns[3:]]

dels = 0
dups = 0

#Get total number of copies of amplicons deleted and duplicated
for ind in vectors.index:
	for repcount in vectors.loc[ind][3:] - normal_cn:
		if repcount > 0:
			dups += repcount
		elif repcount < 0:
			dels -= repcount

print 'Total number of amplicon copies deleted: %i' %(dels)
print 'Total number of amplicon copies duplicated: %i' %(dups)
print

#See how many individuals have only deletions, only duplications, or both
both_inds = 0
del_inds = 0
dup_inds = 0
none_inds = 0

for ind in vectors.index:
	repcounts = vectors.loc[ind][3:] - normal_cn
	if any(x > 0 for x in repcounts) and any(x < 0 for x in repcounts):
		both_inds += 1
	elif any(x > 0 for x in repcounts):
		dup_inds += 1
	elif any(x < 0 for x in repcounts):
		del_inds += 1
	else:
		none_inds += 1

print 'Individuals with no amplicon CNVs: %i out of %i (%i%%)' \
	%(none_inds, len(vectors), float(none_inds)/len(vectors)*100)
print 'Individuals with only amplicon duplications: %i out of %i (%i%%)' \
	%(dup_inds, len(vectors), float(dup_inds)/len(vectors)*100)
print 'Individuals with only amplicon deletions: %i out of %i (%i%%)' \
	%(del_inds, len(vectors), float(del_inds)/len(vectors)*100)
print 'Individuals with both amplicon duplications and deletions: %i out of %i (%i%%)' \
	%(both_inds, len(vectors), float(both_inds)/len(vectors)*100)
print

pie_out_dir = './'

plt.clf()
plt.pie([none_inds, dup_inds, del_inds, both_inds],
		labels=['', '', '', ''],
		colors=['#ffffff', '#D3D3D3', '#696969', '#000000'], 
		wedgeprops = {'linewidth': 0.5, 'edgecolor': 'black'})
plt.axis('equal')
plt.legend(labels=['No CNV', 'Duplication', 'Deletion', 'Complex Variant'], bbox_to_anchor=(1, 1), loc=2)
plt.savefig('%s1000_Genomes_Men_SV_Percent.png' %(pie_out_dir), dpi=300, bbox_inches='tight')
plt.savefig('%s1000_Genomes_Men_SV_Percent.svg' %(pie_out_dir), bbox_inches='tight')

plt.clf()

idx_offset = 1
x_space_used = 0
xtick_locs = []
xtick_labels = []
happcts = {}
#for hapidx, haplogroup in enumerate(sorted(list(set(vectors['Haplogroup'])))):
for hapidx, haplogroup in enumerate(['A', 'B', 'D', 'E', 'C', 'G', 'H', 'I', 'J', 'L', 'T', 'N', 'O', 'Q', 'R']):
	if haplogroup == 'NONE':
		continue
	both_inds = 0
	del_inds = 0
	dup_inds = 0
	none_inds = 0

	hapvectors = vectors.loc[vectors['Haplogroup'] == haplogroup]
	for ind in hapvectors.index:
		repcounts = vectors.loc[ind][3:] - normal_cn
		if any(x > 0 for x in repcounts) and any(x < 0 for x in repcounts):
			both_inds += 1
		elif any(x > 0 for x in repcounts):
			dup_inds += 1
		elif any(x < 0 for x in repcounts):
			del_inds += 1
		else:
			none_inds += 1
	
	print 'Individuals in %s with no amplicon CNVs: %i out of %i (%i%%)' \
		%(haplogroup, none_inds, len(hapvectors), float(none_inds)/len(hapvectors)*100)
	print 'Individuals in %s with only amplicon duplications: %i out of %i (%i%%)' \
		%(haplogroup, dup_inds, len(hapvectors), float(dup_inds)/len(hapvectors)*100)
	print 'Individuals in %s with only amplicon deletions: %i out of %i (%i%%)' \
		%(haplogroup, del_inds, len(hapvectors), float(del_inds)/len(hapvectors)*100)
	print 'Individuals in %s with both amplicon duplications and deletions: %i out of %i (%i%%)' \
		%(haplogroup, both_inds, len(hapvectors), float(both_inds)/len(hapvectors)*100)
	print
	happcts[haplogroup] = float(none_inds)/len(hapvectors)


	plt.pie([none_inds, dup_inds, del_inds, both_inds],
			labels=['', '', '', ''],
			colors=['#ffffff', '#D3D3D3', '#696969', '#000000'],
			wedgeprops={'edgecolor': 'black', 'linewidth': 0.2},
			radius=np.sqrt(len(hapvectors)),
			center=[x_space_used+np.sqrt(len(hapvectors)),0])
	
	
	xtick_locs.append(x_space_used+np.sqrt(len(hapvectors)))
	xtick_labels.append(haplogroup)
	x_space_used += np.sqrt(len(hapvectors)) * 2 + 5

plt.axis('equal')
plt.xticks(xtick_locs, xtick_labels) #FOR DIFFERENT SIZES
plt.tick_params(axis='x', bottom='off')
#plt.subplots_adjust(hspace=.5) #FOR SAME SIZE
plt.legend(labels=['No CNV', 'Duplication', 'Deletion', 'Complex Variant'], loc=2, bbox_to_anchor=(1, 1.5))
plt.savefig('%s1000_Haplogroups_SV_Percent_Scaled.svg' %(pie_out_dir), bbox_inches='tight')
plt.savefig('%s1000_Haplogroups_SV_Percent_Scaled.png' %(pie_out_dir), dpi=300, bbox_inches='tight')















