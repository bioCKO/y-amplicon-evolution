#!/usr/bin/python

'''
Converts Y normalized depth vectors to copy number.
'''

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

import sys
#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/technical_replicates/tech_reps_amplicons.pickle')
#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls.pickle')

#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/gtex/GTEx_Y_amplicon_depths.pickle')

#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls_badctrlfiltered.pickle')

#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons_raw.pickle')

#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons_w100_c11_masked.pickle')



vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_amplicons_calls.pickle')
del vectors['P7']

'''
vectors = vectors.filter([x for x in vectors.columns if x not in complex_and_control])
#vectors = vectors.loc[vectors['Single-copy_depth'] > 1]

for ind in vectors.index:
	for amp in vectors.columns[3:]:
		norm_cov = vectors.loc[ind, amp]
		if np.isnan(norm_cov) or np.isinf(norm_cov):
			vectors.set_value(ind, amp, 0)
		else:
			copy_number = int(norm_cov * copy_numbers[amp] + 0.5)
			vectors.set_value(ind, amp, copy_number)


## vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_calls.pickle')
# import pickle
# outfile = open('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_calls.pickle', 'w')
# pickle.dump(vectors, outfile)
# outfile.close()
# quit()

# import pickle
# outfile = open('/lab/solexa_page/lsteitz/gtex/GTEx_Y_amplicon_calls.pickle', 'w')
# pickle.dump(vectors, outfile)
# outfile.close()
# quit()
# 

import pickle
#outfile = open('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons_raw_calls.pickle', 'w')
outfile = open('/lab/solexa_page/lsteitz/1000_Ys/technical_replicates/tech_reps_amplicons_calls.pickle', 'w')
pickle.dump(vectors, outfile)
outfile.close()
print 'Wrote new output file!'
quit()
'''


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
	
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Depth_Controls/depth_control.svg', bbox_inches='tight')


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

	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Depth_Controls/depth_control_%s.svg' %(hap), bbox_inches='tight')


'''
#Determine relationship between evolutionary distance from reference and number of CNVs

from statsmodels.stats.proportion import proportion_confint
#Gets binomial confidence interval for fraction that have CNVs
from scipy.stats import linregress
from numpy.polynomial.polynomial import polyfit


def cnv_vs_divergence(weight_function, outfile, excluded_amps=[], excluded_haps=[]):
	
	haplogroup_dist_from_ref = {'A': 174.7, 'B': 105.8, 'DE': 76.0, 'C': 75.5, 'G': 54.2, 'H': 54.0, 'IJ': 53.1, 'LTNO': 50.9, 'Q': 35.0, 'R2': 32.9, 'R1': 0.0}
	dist_points = {haplogroup_dist_from_ref[x]:[0,0] for x in haplogroup_dist_from_ref if not x in excluded_haps}


	for ind in vectors.index:
		if vectors.loc[ind, 'Haplogroup'] in ['NONE'] + excluded_haps:
			continue
		if vectors.loc[ind, 'Haplogroup'].startswith('R'):
			for hap in haplogroup_dist_from_ref:
				if vectors.loc[ind, 'Sub-Haplogroup'].startswith(hap):
					dist_points[haplogroup_dist_from_ref[hap]][0] += 1
					break
			for amp in [x for x in vectors.columns[3:] if x not in excluded_amps]:
				if vectors.loc[ind, amp] != copy_numbers[amp]:
					dist_points[haplogroup_dist_from_ref[hap]][1] += 1
					break
		else:
			for hap in haplogroup_dist_from_ref:
				if vectors.loc[ind, 'Haplogroup'] in hap:
					dist_points[haplogroup_dist_from_ref[hap]][0] += 1
					break
			for amp in [x for x in vectors.columns[3:] if x not in excluded_amps]:
				if vectors.loc[ind, amp] != copy_numbers[amp]:
					dist_points[haplogroup_dist_from_ref[hap]][1] += 1
					break

	age_vals = sorted(dist_points.keys())
	cnv_props = [float(dist_points[x][1])/dist_points[x][0] for x in age_vals]
	conf_ints = [proportion_confint(dist_points[x][1], dist_points[x][0], method='beta') for x in age_vals]
	plt.clf()
	plt.errorbar(age_vals, 
					cnv_props, 
					yerr=[[y - x[0] if x[0] > 0 else 0. for x,y in zip(conf_ints,cnv_props)], 
						  [x[1] - y for x,y in zip(conf_ints,cnv_props)]], 
					fmt='s', 
					color='white', 
					ecolor='black')

	cnv_prop_weights = [weight_function(x) for x in conf_ints]
	#cnv_prop_weights = [1] * len(conf_ints)
	#cnv_prop_weights = [1 - (x[1] - max(0,x[0])) for x in conf_ints]
	#cnv_prop_weights = [1/(x[1] - max(0,x[0])) for x in conf_ints]
	#cnv_prop_weights = [1/np.sqrt((x[1] - max(0,x[0]))) for x in conf_ints]
	#print cnv_prop_weights

	#regression_line = linregress(age_vals, cnv_props)
	regression_line_2 = polyfit(age_vals, cnv_props, 1, w=cnv_prop_weights, full=True)

	def calculate_sum_squares(y_data, weights=None):
		y_mean = np.mean(y_data)
		if weights == None:
			weight_vector = [1] * len(y_data)
		else: 
			weight_vector = weights
		return sum([((y - y_mean) * weight)**2 for y, weight in zip(y_data, weight_vector)])

	def calculate_residuals(x_data, y_data, fit_slope, fit_intercept, weights=None):
		if weights == None:
			weight_vector = [1] * len(y_data)
		else: 
			weight_vector = weights
		return sum([((y - ((x  * fit_slope) + fit_intercept)) * weight)**2 for x, y, weight in zip(x_data, y_data, weight_vector)])
		
	def r_squared(x_data, y_data, fit_slope, fit_intercept, weights=None):
		sum_squares = calculate_sum_squares(y_data, weights=weights)
		residuals = calculate_residuals(x_data, y_data, fit_slope, fit_intercept, weights=weights)
		return 1 - (float(residuals)/sum_squares)

	#print regression_line
	#print regression_line_2
	#print r_squared(age_vals, cnv_props, regression_line_2[0][1], regression_line_2[0][0], weights = range(len(cnv_props)))

	r_squared_value = r_squared(age_vals, cnv_props, regression_line_2[0][1], regression_line_2[0][0], weights = cnv_prop_weights)
	
	max_x = max(age_vals) / 10 * 10 + 20
	
	plt.plot([-5, max_x], 
				[-5 * regression_line_2[0][1] + regression_line_2[0][0], 
					max_x * regression_line_2[0][1] + regression_line_2[0][0]],
				'black')
			
	plt.axis([-10,max_x + 10,0,1])
	plt.xlabel('Divergence from Reference (kya)')
	plt.ylabel('Fraction of individuals with CNV')
	plt.text(-5, .94, r'$R^2 = %.3f$' %(r_squared_value))
	#plt.text(-5, .94, r'$R^2 = %.3f$' %(regression_line[2] ** 2))
	#plt.text(-5, .90, r'$p = %f$' %(regression_line[3]))
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/%s.png' %(outfile), dpi=200)

cnv_vs_divergence(lambda x: 1, 'age_vs_percent_CNV_unweighted')
cnv_vs_divergence(lambda x: 1/np.sqrt((x[1] - max(0,x[0]))), 'age_vs_percent_CNV_weighted')
cnv_vs_divergence(lambda x: 1, 'age_vs_percent_CNV_no_AB_unweighted', excluded_haps=['A','B'])
cnv_vs_divergence(lambda x: 1/np.sqrt((x[1] - max(0,x[0]))), 'age_vs_percent_CNV_no_AB__weighted', excluded_haps=['A','B'])
cnv_vs_divergence(lambda x: 1, 'age_vs_percent_CNV_no_ABDN_unweighted', excluded_haps=['A','B', 'D', 'N'])
cnv_vs_divergence(lambda x: 1/np.sqrt((x[1] - max(0,x[0]))), 'age_vs_percent_CNV_no_ABDN_weighted', excluded_haps=['A','B', 'D', 'N'])
'''

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

pie_out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Pie Charts/'
'''
plt.clf()
plt.pie([none_inds, dup_inds, del_inds, both_inds],
		labels=['', '', '', ''],
		colors=['#ffffff', '#D3D3D3', '#696969', '#000000'], 
		wedgeprops = {'linewidth': 0.5, 'edgecolor': 'black'})
plt.axis('equal')
plt.legend(labels=['No CNV', 'Duplication', 'Deletion', 'Complex Variant'], bbox_to_anchor=(1, 1), loc=2)
plt.savefig('%s1000_Genomes_Men_SV_Percent.png' %(pie_out_dir), dpi=300, bbox_inches='tight')
plt.savefig('%s1000_Genomes_Men_SV_Percent.svg' %(pie_out_dir), bbox_inches='tight')
'''
plt.clf()
'''
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
'''
# 	if haplogroup == 'K':
# 		idx_offset = 0
# 		continue
# 	plt.subplot(4, 4, hapidx + idx_offset)
# 	plt.pie([none_inds, dup_inds, del_inds, both_inds],
# 			labels=['', '', '', ''],
# 			colors=['#ffffff', '#D3D3D3', '#696969', '#000000'],
# 			wedgeprops={'linewidth': 0.2},)
# 	
# 	plt.axis('equal')
# 	plt.xlabel(haplogroup)
'''

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
#plt.savefig('%s1000_Haplogroups_SV_Percent_Scaled.svg' %(pie_out_dir), bbox_inches='tight')
#plt.savefig('%s1000_Haplogroups_SV_Percent_Scaled.png' %(pie_out_dir), dpi=300, bbox_inches='tight')
'''














