#!/usr/bin/python

'''
Simulates mutation over the Y-chromosome phylogenetic tree and performs the branch sorting 
analysis on each simulation.
'''

import sys
import ete3
import copy
import pickle
import random
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from scipy.stats import kstest
from branch_sorting import fitch_up, fitch_down


def evolve_tree(newick, mut_rate=0):
	t = copy.deepcopy(newick)

	gen_per_snp = 1 / (25 * .0076) #Generations per unit of branch length/SNP; .0076 is SNPs per year, 25 is years per generation
	for node in t.traverse():
		node.add_feature('mutation', 'no')

		rand_draw = random.random()
		if rand_draw >= (1 - mut_rate) ** (node.dist * gen_per_snp):
			node.mutation = 'yes'
		
	return t


parser = argparse.ArgumentParser(description = "Calculates a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree")

parser.add_argument("-t", "--tree", help = "Phylogenetic tree of Y chromosomes in newick format")
parser.add_argument("-o", "--output_directory", default = "./", help = "Directory to write output files to (default: current directory)")
parser.add_argument("-m", "--mutation_rate", default = .000383, help = "Amplicon CNV mutation rate (default = 3.83e-4)")
parser.add_argument("-v", "--variant_annotation", help = "File containing individuals annotated with amplicon variants")
parser.add_argument("-s", "--simulated_pvals", action = "store_true", help = "Create file listing p-values of 1000 simulations")

args = parser.parse_args()

tree = ete3.Tree(args.tree, format=1)


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

	
plt.clf()
plt.figure(figsize=[7,7])
plt.xlabel('Evolutionary time\n(SNPs)')
plt.ylabel('Fraction of variants observed')
ks_ps = []
mut_events = []


for sim in xrange(10):
	sim_tree = evolve_tree(tree, mut_rate=args.mutation_rate)
	
	branches = []
	for branch in sim_tree.traverse():
	
		distance = 0
		nextdist = branch.dist
		mutation = branch.mutation
	
		if branch.is_leaf():
			leaf_status = 1
		else:
			leaf_status = 0
	
		#distance from leaves
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
		
		branches.append([distance, nextdist, mutation, leaf_status])


	branches.sort(key=lambda x: x[0], reverse=True)
	
	mut_events.append(len([x[2] for x in branches if x[2].startswith('yes')]))

	age = 0
	mut_times = []
	for branch in branches:
		#Only a single mutation can occur per branch
		if branch[2] == 'yes':
			mut_times.append(age + branch[1]/2.)
		
		age += branch[1]
		
	
	if len(mut_times) > 0:
		ks_ps.append(kstest(np.array(mut_times)/float(age), 'uniform'))
	else:
		ks_ps.append((0,1))
	

	if len(mut_times) > 0:
		plt.plot(mut_times, np.arange(1, len(mut_times) + 1)/float(len(mut_times)), 'gray', alpha=0.3)


plt.plot([0,age],[0,1], 'r', lw=1.5)
plt.axis([0, age, 0, 1])


plt.savefig('%ssimulated_branch_sorting.png' %(args.output_directory), dpi=300)
plt.savefig('%ssimulated_branch_sorting.svg' %(args.output_directory))

if args.simulated_pvals:
	import pickle
	outfile = open('%ssimulated_branch_sorting_pvalues.pickle' %(args.output_directory), 'w')
	pickle.dump(ks_ps, outfile)
	outfile.close()