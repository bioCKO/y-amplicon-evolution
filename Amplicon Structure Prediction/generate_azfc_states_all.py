#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Amplicon:
	def __init__(self, name, orientation):
		self.name = name
		if orientation not in ['f', 'r']:
			raise TypeError('Orientation must be \'f\' or \'r\'')
		self.orientation = orientation
		
	def invert(self):
		if self.orientation == 'r':
			return Amplicon(self.name, 'f')
		elif self.orientation == 'f':
			return Amplicon(self.name, 'r')
			
	def __str__(self):
		return str(self.__dict__)

	def __eq__(self, other): 
		return self.__dict__ == other.__dict__
		
		
def delete_amplicon(azfc, ind_1, ind_2):
	if azfc[ind_1].name == azfc[ind_2].name and azfc[ind_1].orientation == azfc[ind_2].orientation:
		return azfc[:ind_1] + azfc[ind_2:]
	return azfc

def duplicate_amplicon(azfc, ind_1, ind_2):
	if azfc[ind_1].name == azfc[ind_2].name and azfc[ind_1].orientation == azfc[ind_2].orientation:
		return azfc[:ind_2] + azfc[ind_1:ind_2] + azfc[ind_2:]
	return azfc
		
def invert_amplicon(azfc, ind_1, ind_2):
	if azfc[ind_1].name == azfc[ind_2].name and azfc[ind_1].orientation != azfc[ind_2].orientation:
		return azfc[:ind_1 + 1] + [x.invert() for x in azfc[ind_1 + 1:ind_2][::-1]] + azfc[ind_2:]
	return azfc
	
def count_amplicons(azfc):
	amps = ['IR1', 'IR5', 'blue', 'teal', 'green', 'red', 'gray', 'yellow']
	#This order is important, since the counts will be returned as a tuple with 
	#no other identification
	amp_counts = [len([y for y in azfc if y.name == x]) for x in amps]
	amp_counts[0] += 1 #Add the IR1 copy not in AZFc
	amp_counts[1] += 2 #Add the IR5 copies not in AZFc
	return tuple(amp_counts)
	
def plot_azfc(azfc, outfile):
	amp_lens = {'IR1':0.2, 'IR5':0, 'blue':3, 'teal':2, 'green':5, 'red':2, 'gray':1.5, 'yellow':10}
	loc = 0
	head_len = 1
	body_width = 0.1
	head_width = 0.15
	plt.clf()
	for amp in azfc:
		if amp_lens[amp.name] > 0:
			if amp.name.startswith('IR'):
				plt.arrow(loc, 0, amp_lens[amp.name], 0,
							width=body_width, head_width=0, head_length=0, fc='white')
			elif amp.orientation == 'f':
				plt.arrow(loc, 0, amp_lens[amp.name] - head_len, 0, 
							width=body_width, head_width=head_width, head_length=head_len, fc=amp.name)
			elif amp.orientation == 'r':
				plt.arrow(loc + amp_lens[amp.name], 0, -(amp_lens[amp.name] - head_len), 0, 
							width=body_width, head_width=head_width, head_length=head_len, fc=amp.name)
		loc += amp_lens[amp.name]
	plt.axis([-1, loc + 1, -1, 1])
	plt.axis('off')
	plt.savefig(outfile, dpi=200)
	
def plot_multiple_azfc(azfclist, outfile):
	amp_lens_props = {'IR1':0.2, 'IR5':0, 'blue':3, 'teal':2, 'green':5, 'red':2, 'gray':1.5, 'yellow':10}
	amp_lens_modifier = 3
	amp_lens = {x: amp_lens_props[x] * amp_lens_modifier for x in amp_lens_props}
	head_len = 3
	body_width = 0.3
	head_width = 0.45
	max_len = 0
	plt.clf()
	for azfc_idx, azfc in enumerate(azfclist):
		loc = 0
		for amp in azfc:
			if amp_lens[amp.name] > 0:
				if amp.name.startswith('IR'):
					plt.arrow(loc, -azfc_idx, amp_lens[amp.name], 0,
								width=body_width, head_width=0, head_length=0, fc='white', lw=0.5)
				elif amp.orientation == 'f':
					plt.arrow(loc, -azfc_idx, amp_lens[amp.name] - head_len, 0, 
								width=body_width, head_width=head_width, head_length=head_len, fc=amp.name, lw=0.5)
				elif amp.orientation == 'r':
					plt.arrow(loc + amp_lens[amp.name], -azfc_idx, -(amp_lens[amp.name] - head_len), 0, 
								width=body_width, head_width=head_width, head_length=head_len, fc=amp.name, lw=0.5)
			loc += amp_lens[amp.name]
		if loc > max_len:
			max_len = loc
	plt.axis([-1, max_len + 1, -azfc_idx -1, 1])
	plt.axis('off')
	plt.savefig(outfile, dpi=500)

def find_closest_predicted(copy_number, cn_predict):
	comps = []
	for test_copy_number in cn_predict:
		if test_copy_number == 'other':
			continue
		comps.append([test_copy_number, sum([abs(x - y) for x, y in zip(copy_number, test_copy_number)])])
	comps.sort(key=lambda x: x[1])
	return comps[:5]

def least_necessary_rearrangements(copy_number, cn_predict):
	if copy_number in cn_predict:
		if len(cn_predict[copy_number]['reference']) > 0:
			return 0
		elif len(cn_predict[copy_number]['one_nahr']) > 0:
			return 1
		elif len(cn_predict[copy_number]['two_nahr']) > 0:
			return 2
		elif len(cn_predict[copy_number]['three_nahr']) > 0:
			return 3
	else:
		return 4

def least_necessary_rearrangements_structure(copy_number, cn_predict):
	if copy_number in cn_predict:
		if len(cn_predict[copy_number]['reference']) > 0:
			return cn_predict[copy_number]['reference']
		elif len(cn_predict[copy_number]['one_nahr']) > 0:
			return cn_predict[copy_number]['one_nahr']
		elif len(cn_predict[copy_number]['two_nahr']) > 0:
			return cn_predict[copy_number]['two_nahr']
		elif len(cn_predict[copy_number]['three_nahr']) > 0:
			return cn_predict[copy_number]['three_nahr']
	else:
		return None
		
if __name__ == '__main__':
	
	named_structures = {(2, 4, 4, 2, 3, 4, 2, 2):'reference',
						(2, 3, 3, 2, 2, 2, 1, 1):'gr/gr deletion',
						(1, 3, 3, 2, 1, 2, 2, 1):'b2/b3 deletion',
						(2, 5, 5, 2, 4, 6, 3, 3):'gr/gr duplication',
						(3, 5, 5, 2, 5, 6, 2, 3):'b2/b3 duplication'}
	
	import itertools
	from collections import defaultdict
	import pandas as pd
	
	amp_names = ['blue', 'green', 'red', 'gray', 'yellow'] #Teal and IR1 not included since they can't mediate NAHR
	
	#Reference orientation
	azfc_reference = [Amplicon('blue', 'r'),
					  Amplicon('teal', 'r'),
					  Amplicon('teal', 'f'),
					  Amplicon('blue', 'f'),
					  Amplicon('IR1', 'f'),
					  Amplicon('green', 'r'),
					  Amplicon('red', 'r'),
					  Amplicon('red', 'f'),
					  Amplicon('gray', 'r'),
					  Amplicon('blue', 'r'),
					  Amplicon('yellow', 'r'),
					  Amplicon('IR5', 'r'),
					  Amplicon('green', 'r'),
					  Amplicon('red', 'r'),
					  Amplicon('red', 'f'),
					  Amplicon('green', 'f'),
					  Amplicon('IR5', 'f'),
					  Amplicon('yellow', 'f'),
					  Amplicon('blue', 'f'),
					  Amplicon('gray', 'f')]


	#Orientations possible through 1 NAHR event
	
	one_nahr_muts = []
	for amp in amp_names:
		amp_indices = [x for x in xrange(len(azfc_reference)) if azfc_reference[x].name == amp]
		for amp_pair in itertools.combinations(amp_indices, 2):
			deletion = delete_amplicon(azfc_reference, *amp_pair)
			duplication = duplicate_amplicon(azfc_reference, *amp_pair)
			inversion = invert_amplicon(azfc_reference, *amp_pair)
			for variant in [deletion, duplication, inversion]:
				if variant not in [azfc_reference] + one_nahr_muts:
					one_nahr_muts.append(variant)
	
	

	#Orientations possible through 2 NAHR events
	
	two_nahr_muts = []
	for azfc in one_nahr_muts:
		for amp in amp_names:
			amp_indices = [x for x in xrange(len(azfc)) if azfc[x].name == amp]
			for amp_pair in itertools.combinations(amp_indices, 2):
				deletion = delete_amplicon(azfc, *amp_pair)
				duplication = duplicate_amplicon(azfc, *amp_pair)
				inversion = invert_amplicon(azfc, *amp_pair)
				for variant in [deletion, duplication, inversion]:
					if variant not in [azfc_reference] + one_nahr_muts + two_nahr_muts:
						two_nahr_muts.append(variant)
	
	
	#Orientations possible through 3 NAHR events
	
	three_nahr_muts = []
	for azfc in two_nahr_muts:
		for amp in amp_names:
			amp_indices = [x for x in xrange(len(azfc)) if azfc[x].name == amp]
			for amp_pair in itertools.combinations(amp_indices, 2):
				deletion = delete_amplicon(azfc, *amp_pair)
				duplication = duplicate_amplicon(azfc, *amp_pair)
				inversion = invert_amplicon(azfc, *amp_pair)
				for variant in [deletion, duplication, inversion]:
					if variant not in [azfc_reference] + one_nahr_muts + two_nahr_muts + three_nahr_muts:
						three_nahr_muts.append(variant)
	
	#Orientations possible through non-amplicon-NAHR events
	
	non_nahr_muts = []
	for dstart in xrange(0, len(azfc_reference)):
		for dend in xrange(dstart + 1, len(azfc_reference) + 1):
			non_nahr_muts.append(azfc_reference[:dstart] + azfc_reference[dend:])
			non_nahr_muts.append(azfc_reference[:dend] + azfc_reference[dstart:])

	#Make dict of copy number states of non-NAHR events
	
	non_nahr_cns = {}
	for mut in non_nahr_muts:
		if not count_amplicons(mut) in non_nahr_cns:
			non_nahr_cns[count_amplicons(mut)] = []
		non_nahr_cns[count_amplicons(mut)].append(mut)

	
	#Make dict of copy number states and the AZFc architectures with those copy numbers
	
	copy_numbers = defaultdict(lambda: {'reference': [], 'one_nahr': [], 'two_nahr': [], 'three_nahr': []})
	
	copy_numbers[count_amplicons(azfc_reference)]['reference'].append(azfc_reference)
	
	for mut in one_nahr_muts:
		copy_numbers[count_amplicons(mut)]['one_nahr'].append(mut)
	
	for mut in two_nahr_muts:
		copy_numbers[count_amplicons(mut)]['two_nahr'].append(mut)
	
	for mut in three_nahr_muts:
		copy_numbers[count_amplicons(mut)]['three_nahr'].append(mut)
	
	print 'Predicted number of AZFc CN states: %i' %(len(copy_numbers))
	print len(one_nahr_muts), len(two_nahr_muts), len(three_nahr_muts)
	print 'Predicted number of AZFc structures: %i' %(sum([sum([len(x) for x in y.values()]) for y in copy_numbers.values()]))
	print
	#Find copy number states from 1000 Genomes data
	
	
	
	out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Observed Variants/'
	
	'''
	present_states = [[x, len(azfc_cn_states[x])] for x in azfc_cn_states if len(azfc_cn_states[x]) > 0]
	present_states.sort(key=lambda x: x[1], reverse=True)
	
	present_states_structures = []
	
	for state in present_states:
		if state[0] == 'other':
			continue
		for structure in least_necessary_rearrangements_structure(state[0], copy_numbers):
			present_states_structures.append(structure)
		present_states_structures.append([])
	'''
	
	#plot_multiple_azfc(two_nahr_muts, '%sTwo Step States.svg' %(out_dir))
	
	#quit()
	
	
	azfc_cn_states = {x:[] for x in copy_numbers}
	azfc_cn_states['other'] = {}
	
	import sys
	vectors = pd.read_pickle(sys.argv[1])
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_calls.pickle')
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_calls_badctrlfiltered.pickle')
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons_properly_masked_calls.pickle')
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/properly_masked_vectors/1000_Genomes_amplicons_w100_c11_masked_calls.pickle')
	
	#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/gtex/GTEx_Y_amplicon_calls.pickle')
	
	vectors = vectors.drop(['IR3', 'P8', 'P7', 'P6', 'P5', 'P4', 'IR2'], 1)
	
	
	others_haplogroups = defaultdict(lambda: [])
	others_men = defaultdict(lambda: [])
	for ind in vectors.index:
		ind_state = tuple(map(int,vectors.loc[ind][3:]))
		if ind_state in azfc_cn_states:
			azfc_cn_states[ind_state].append(ind)
		else:
			azfc_cn_states['other'][ind] = ind_state
			others_haplogroups[ind_state].append(vectors.loc[ind]['Sub-Haplogroup'])
			others_men[ind_state].append(ind)
	
	
	azfc_cn_haplogroups = {x:[] for x in copy_numbers}
	
	
	number_of_rearrangements = [0,0,0,0,0]
	for ind in vectors.index:
		ind_state = tuple(map(int,vectors.loc[ind][3:]))
		ind_rearrangements = least_necessary_rearrangements(ind_state, copy_numbers)
		number_of_rearrangements[ind_rearrangements] += 1
	
	print 'Total men: %i' %(sum(number_of_rearrangements))
	for recount in xrange(5):
		print 'Men with %i rearrangements: %i' %(recount, number_of_rearrangements[recount])
	print
	
	others = {}
	for thing in azfc_cn_states['other'].values():
		if thing not in others:
			others[thing] = 1
		else:
			others[thing] += 1        	
	others = [[a, others[a]] for a in others]
	others.sort(key=lambda x: x[1], reverse=True)
	
	found_states = []
	for x in azfc_cn_states:
		if len(azfc_cn_states[x]) > 0:
			found_states.append([x, len(azfc_cn_states[x]), least_necessary_rearrangements(x, copy_numbers)])
	print 'CN states in 1000 Genomes men:'
	for x in sorted(found_states, key=lambda x: x[1], reverse=True):
		print x[0],
		if x[0] in named_structures:
			print '(%s)' %(named_structures[x[0]]),
		print '\t# of men: %i\tMin # of rearrangements: %i' %(x[1], x[2]),
		if x[0] in non_nahr_cns and x[2] > 1:
			print '(or non-NAHR)'
		else:
			print
	print
	'''	
	amps = ['IR1', 'IR5', 'blue', 'teal', 'green', 'red', 'gray', 'yellow']
	simp_indels = 0
#	print 'Unpredicted states:'
	for x in others:
#		print '%s\t# of men: %i' %(str(x[0]), x[1]),

#		print ' (%s)' %(', '.join(['%s: %i' %(h, others_haplogroups[x[0]].count(h)) for h in sorted(list(set(others_haplogroups[x[0]])))])),

		non_nahr_off = find_closest_predicted(x[0], non_nahr_cns)[0][1]
#		print '\t(%i of off non-NAHR)' %(non_nahr_off)
		if non_nahr_off == 0:
			simp_indels += 1
#		print '\tClosest predicted structures:'
		for prestruct in find_closest_predicted(x[0], copy_numbers)[:2]:
#			print '\t%s%s\t' %(str(prestruct[0]), ' ' + named_structures[prestruct[0]] if prestruct[0] in named_structures else ''),
			diff_str = []
			for amp_idx in xrange(len(x[0])):
				if x[0][amp_idx] != prestruct[0][amp_idx]:
					diff = x[0][amp_idx] - prestruct[0][amp_idx]
					diff_str.append('%s%i %s' %('+' if diff > 0 else '-', abs(diff), amps[amp_idx]))			
#			print 'CN difference: %i (%s)' %(prestruct[1], ', '.join(diff_str))	
#		print		
#	print
	'''


	amps = ['IR1', 'IR5', 'blue', 'teal', 'green', 'red', 'gray', 'yellow']
	simp_indels = 0
	print 'Unpredicted states:'
	for x in others:
		print '%s\t# of men: %i' %(str(x[0]), x[1]),

		print ' (%s)' %(', '.join(['%s: %i' %(h, others_haplogroups[x[0]].count(h)) for h in sorted(list(set(others_haplogroups[x[0]])))])),
		print ' (%s)' %(', '.join(others_men[x[0]])),

		non_nahr_off = find_closest_predicted(x[0], non_nahr_cns)[0][1]
		print '\t(%i of off non-NAHR)' %(non_nahr_off)
		if non_nahr_off == 0:
			simp_indels += 1
		print '\tClosest predicted structures:'
		for prestruct in find_closest_predicted(x[0], copy_numbers)[:2]:
			print '\t%s%s\t' %(str(prestruct[0]), ' ' + named_structures[prestruct[0]] if prestruct[0] in named_structures else ''),
			diff_str = []
			for amp_idx in xrange(len(x[0])):
				if x[0][amp_idx] != prestruct[0][amp_idx]:
					diff = x[0][amp_idx] - prestruct[0][amp_idx]
					diff_str.append('%s%i %s' %('+' if diff > 0 else '-', abs(diff), amps[amp_idx]))			
			print 'CN difference: %i (%s)' %(prestruct[1], ', '.join(diff_str))	
		print		
	print
	





	print 'Others explainable by a non-NAHR event: %i out of %i' %(simp_indels, len(others))
	off_by_one = [x for x in others if find_closest_predicted(x[0], copy_numbers)[0][1] == 1]
	print 'Others off by one: %i out of %i' %(len(off_by_one), len(others))
	
	
	#Plot all of the observed CNVs
	
	out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Observed Variants/'
	
	present_states = [[x, len(azfc_cn_states[x])] for x in azfc_cn_states if len(azfc_cn_states[x]) > 0]
	present_states.sort(key=lambda x: x[1], reverse=True)
	
	present_states_structures = []
	
	for state in present_states:
		if state[0] == 'other':
			continue
			
		#print azfc_cn_states[state[0]]
			
		for structure in least_necessary_rearrangements_structure(state[0], copy_numbers):
			present_states_structures.append(structure)
		present_states_structures.append([])
	
	#plot_multiple_azfc(present_states_structures, '%sDetected States.svg' %(out_dir))
	
	
	
	
	
	
	
	
	
	
	
	