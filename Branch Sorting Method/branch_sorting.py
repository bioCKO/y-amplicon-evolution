#!/usr/bin/python

'''
Calculates a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree.
Does so by sorting all branches of the tree from oldest to youngest, and then compares 
the CDF of mutation events over time to a null model where those events are expected to 
occur neutrally (that is, same probability at any time throughout the evolutionary history).
Using a KS test between the observed and null CDFs, we get a p-value.
'''

import ete3
import pickle
import random
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from scipy.stats import kstest


def fitch_up(tree): #Does first stage of Fitch algorithm, going up the tree
	event_count = 0

	for node in tree.traverse('postorder'):
		if node.is_leaf():
			continue
		states_1, states_2 = [x.cnv_state for x in node.children]
		if not type(states_1) is set:
			states_1 = set([states_1])
		if not type(states_2) is set:
			states_2 = set([states_2])
		if states_1 & states_2:
			node.add_feature('cnv_state', states_1 & states_2)
		else:
			node.add_feature('cnv_state', states_1 | states_2)
			event_count += 1
	
	return event_count

def fitch_down(tree): #Does second stage of Fitch algorithm, going down the tree
	for node in tree.traverse('preorder'):
		if node.is_leaf() or node.is_root():
			continue
		if node.up.cnv_state < node.cnv_state:
			node.cnv_state = node.up.cnv_state


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Calculates a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree")
	
	parser.add_argument("-t", "--tree", help = "Phylogenetic tree of Y chromosomes in newick format")
	parser.add_argument("-o", "--output_directory", default = "./", help = "Directory to write output files to (default: current directory)")
	parser.add_argument("-v", "--variant_annotation", help = "File containing individuals annotated with amplicon variants")
	parser.add_argument("-n", "--no_plots", action = "store_true", help = "Do not create plots of real and shuffled branches")
	parser.add_argument("-s", "--shuffled_pvals", action = "store_true", help = "Create file listing p-values of 1000 shuffles")

	args = parser.parse_args()

	tree = ete3.Tree(args.tree, format=1)

	#Annotate tree with variants
	yvarfile = open(args.variant_annotation, 'r')
	yvars = {}
	for line in yvarfile:
		if not line.startswith('#'):
			data = line.split('\t')		
			yvars[data[0]] = data[1]


	good_leaves = []
	for leaf in tree.iter_leaves():
		leafname = leaf.name.split('_')[0]
		if leafname in yvars:
			good_leaves.append(leaf.name)
	tree.prune(good_leaves)

	#Men with a CNV of a single copy of a single amplicon without evidence of a CNV from partial CNV or FISH analyses, used to analyze "high-confidence" dataset
	one_amp_list = ['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA19774', 'NA18982', 'NA19198', 'NA19451', 'HG04219', 'HG03021', 'HG03015', 'NA18977']

	#Men with a CNV of a single copy of a single amplicon, including even those with evidence of a CNV from partial CNV or FISH analyses. Not shown in paper, but this "super-high-confidence" dataset also returns results that differ only minimally from what is shown
	#one_amp_list['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG02512', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA21093', 'HG01842', 'NA21112', 'HG01767', 'HG01530', 'NA21111', 'NA19774', 'NA18982', 'NA19198', 'NA12748', 'HG03660', 'NA19451', 'HG04219', 'NA20884', 'HG03021', 'HG03015', 'NA18977', 'NA20536']


	for leaf in tree.iter_leaves():
		leafname = leaf.name.split('_')[0]
		if leafname in yvars:
			if yvars[leafname] == 'v1': #Use this line to count all amplicon CNVs
			#if yvars[leafname] == 'v1' or leafname in one_amp_list: #Use this line to exclude men with a CNV of a single copy of a single amplicon
				leaf.add_feature('cnv_state', set(['reference']))
			else:
				leaf.add_feature('cnv_state', set([yvars[leafname]]))

	#Name unnamed nodes for ease of annotation later
	branch_no = 1
	for branch in tree.traverse():
		if branch.name == '':
			branch.name = str(branch_no)
			branch_no += 1

	#Use Fitch's algorithm to find locations of mutation events
	fitch_events = fitch_up(tree)
	fitch_down(tree)


	# #This commented-out section is used to make figures of locations in the tree where the 
	# #location of a mutation event is ambiguous based on Fitch's algorithm. The list of
	# #manually curated branches in the section below is based on visual inspection of the 
	# #figures generated with this section.
	# no_states = []
	# ctcount = 0
	# for branch in tree.traverse():
	# 	no_states.append([branch.cnv_state, len(branch.cnv_state)])
	# 	if len(branch.cnv_state) > 1:
	# 		for node in branch.traverse():
	# 			if node.is_leaf():
	# 				nodename = node.name.split('_')[0]
	# 				node.add_face(ete3.TextFace(' %s' %(yvars[nodename])), 1, 'branch-right')
	# 			else:
	# 				node.add_face(ete3.TextFace(node.name), 1, 'branch-right')
	# 
	# 	
	# 		bl_style = ete3.TreeStyle()
	# 		bl_style.show_branch_length = True
	# 		ctcount += 1
	# 		branch.render('%scomplex_trees_%i.pdf' %(args.output_directory, ctcount), tree_style=bl_style)


	#Annotate branches in which mutation occurred:
	for branch in tree.traverse():
	
		if branch.is_root():
			branch.add_feature('mutation', 'no')
			continue
	
	
		#manually set mutation states for one subtree where two adjacent gr/gr deletions make it think it's going from del to reference
		if branch.name in ['106|93', '100|3']:
			branch.add_feature('mutation', 'no')
		elif branch.name in ['HG00442_O2', '105|53']:
			branch.add_feature('mutation', 'yes')
	
		#This is the list of manually curated complex parts of the tree, 
		#based on my knowledge of the actual variant structures.
		elif branch.name in ['HG03812_H', 
							 '20|94',
							 '266|N-M231|175',
							 '257|77',
							 '264|26',
							 'HG02696_Q1a',
							 '143|8',
							 'HG01097_G2', #remove for one amp filtering
							 'NA12546_G2',
							 'HG03615_J2',
							 '91|55',
							 'HG03908_J2',
							 '41|23',
							 'HG01982_R1a',
							 'HG01617_R1a',
							 '137',
							 'HG03745_H1',
							 'HG02070_O3',
							 'NA18982_O3', #remove for one amp filtering
							 'HG02410_O2',
							 'HG02408_O2',
							 'HG02394_O2',
							 'HG03196_E1b',
							 'NA18504_E1b']:
	
			branch.add_feature('mutation', 'yes')
	
		elif branch.name in ['246|16'] and branch.cnv_state != set(['reference']): #There are two branches with this name; only annotate the one with a mutation
			branch.add_feature('mutation', 'yes')
	
		#Unambiguous mutations
		elif len(branch.cnv_state) == 1 and len((branch.up).cnv_state) == 1 and branch.cnv_state != (branch.up).cnv_state:
			branch.add_feature('mutation', 'yes')
	
		else:
			branch.add_feature('mutation', 'no')

		
	#Calculate each branch's age
	branches = []
	for branch in tree.traverse():
	
		distance = 0
		nextdist = branch.dist
		mutation = branch.mutation
	
		if branch.is_leaf():
			leaf_status = 1
		else:
			leaf_status = 0
	
	
		#average distance from leaves
		leaf_dists = []
		for leaf in branch.iter_leaves():
			leaf_dist = branch.dist/2. #This metric needs to include the node's length itself, because I'm really looking at branches, not nodes.
			#The division by 2 is because, for longer branches, I want to consider the "average" time they represent.
			#For example, If I count from the bottom, that treats a leaf with a long way to the preceding split 
			#and a short way the same, which seems wrong. On the other hand, if I count from the top, that treats a
			#leaf with a long way the same as a much higher branch that has the same distance from the leaves the same, 
			#which also seems wrong; a mutation in the former could have taken place over a long period of time, from very 
			#recently to long ago; in the latter, it must have happened long ago. Dividing by 2 is a simple way to account for this.
			node = leaf
			while not node == branch:
				leaf_dist += node.dist
				node = node.up
		
			leaf_dists.append(leaf_dist)
		distance = np.mean(leaf_dists)
	
	
		branches.append([distance, nextdist, mutation, leaf_status, branch.name])


	#Sort branches by age
	branches.sort(key=lambda x: (x[0], x[4]), reverse=True)

	#Calculate how many mutation events are in the first half of the tree and in the entire tree
	tot_age = sum([x[1] for x in branches])
	first_half_muts = 0
	age = 0
	branch_idx = 0
	while age < tot_age/2:
		if branches[branch_idx][2] == 'yes':
			first_half_muts += 1
		age += branches[branch_idx][1]
		branch_idx += 1

	print 'Mutation events in first half of tree: %i' %(first_half_muts)

	while age < tot_age:
		if branches[branch_idx][2] == 'yes':
			first_half_muts += 1
		age += branches[branch_idx][1]
		branch_idx += 1

	print 'Mutation events in entire tree: %i' %(first_half_muts)


	#Make plots of sorted and shuffled branches and calculate KS p-value
	if not args.no_plots:
		plt.clf()

		plt.figure(figsize=[7,7])
		age = 0
		zmut_times = []
		for branch in branches:
			if branch[2] == 'yes':
				zmut_times.append(age + branch[1]/2.)
			age += branch[1]

		ks_p = kstest(np.array(zmut_times)/float(age), 'uniform')
		print 'KS test p-value: %.3e' %(ks_p[1])
		print 'Total branch length of tree: %.1f' %(age)


	
		#Make figure real data and 1000 shuffles
		plt.clf()
		ks_ps = []
		for subfig in xrange(1,1001):
			random.shuffle(branches)
			age = 0
			mut_times = []
			for branch in branches:
				if branch[2] == 'yes':
					mut_times.append(age + branch[1]/2.)
				age += branch[1]

			ks_ps.append(kstest(np.array(mut_times)/float(age), 'uniform'))
			plt.plot(mut_times, np.arange(0,1,1./len(mut_times)), 'gray', alpha=0.3)


		branches.sort(key=lambda x: (x[0], x[4]), reverse=True)
		age = 0
		mut_times = []
		for branch in branches:
			if branch[2] == 'yes':
				mut_times.append(age + branch[1]/2.)
			age += branch[1]


		plt.plot(mut_times, np.arange(0,1,1./len(mut_times)), 'b', lw=1.5)
		plt.plot([0,age],[0,1], 'r', lw=1.5)
		plt.axis([0, age, 0, 1])
		plt.xlabel('Evolutionary time\n(SNPs)')
		plt.ylabel('Fraction of variants observed')
		plt.title(r'$p = %s\times10^{-%i}$ (KS test)' %(('%3e' %ks_p[1])[:5], int(str(ks_p[1]).split('-')[-1])))

		out_name = 'branch_sorting'
		plt.savefig('%s%s.png' %(args.output_directory, out_name), dpi=300)
		plt.savefig('%s%s.svg' %(args.output_directory, out_name))

	#Save pickled list of KS test results for 1000 shuffles
	if args.shuffled_pvals:
		outfile = open('%s/shuffled_ks_pvals.pickle' %(args.output_directory), 'w')
		pickle.dump(ks_ps, outfile)
		outfile.close()
