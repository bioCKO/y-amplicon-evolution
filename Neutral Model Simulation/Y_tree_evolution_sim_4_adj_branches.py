#!/nfs/apps/anaconda2/bin/python

'''
Hashbang for running on local machine: #!/usr/local/bin/miniconda2/bin/python
'''
'''
Simulates amplicon mutation in 1000 Ys tree to get haplogroup CNV % distributions
'''

import ete3
import random
import copy


def evolve_node(node, dup_mut=0.00043, del_mut=0.0, dup_rev=0.00086, del_rev=0.0):
	#dup_mut = Duplication mutation rate (events per generation)
	#del_mut = Deletion mutation rate (events per generation)
	#dup_rev = Duplication reversion rate (events per generation)
	#del_rev = Deletion reversion rate (events per generation)
	
	rand_draw = random.random()
	if node.cnv_state == 'reference':
		if rand_draw < dup_mut:
			node.cnv_state = 'duplication'
			return
		elif rand_draw < dup_mut + del_mut:
			node.cnv_state = 'deletion'
			return
	elif node.cnv_state == 'duplication':
		if rand_draw < dup_rev:
			node.cnv_state = 'reference'
			return
	elif node.cnv_state == 'deletion':
		if rand_draw < del_rev:
			node.cnv_state = 'reference'
			return
	else:
		print 'ERROR: UNRECOGNIZED CNV STATE'
		quit()


def evolve_tree(newick, dup_mut=0.00043, del_mut=0.0, dup_rev=0.00086, del_rev=0.0):
	t = copy.deepcopy(newick)
	
	mutstyle = ete3.NodeStyle()
	mutstyle['fgcolor'] = 'red'

	gen_per_snp = 1 / (25 * .0076) #Generations per unit of branch length/SNP; .0076 is SNPs per year, 25 is years per generation
	for node in t.traverse():
		if node.is_root():
			node.add_feature('cnv_state', 'reference')
			continue
		node.add_feature('cnv_state', node.up.cnv_state)
		orig_state = node.cnv_state
		
		###
		##Full generational simulation; use when multiple mutations per branch are allowed
		#for generation in xrange(int(node.dist * gen_per_snp)):
		#	evolve_node(node, dup_mut=dup_mut, del_mut=del_mut, dup_rev=dup_rev, del_rev=del_rev)
		###
		
		###
		#Shortcut simulation; more accurate times since generations aren't rounded down
		#Right now, does NOT WORK for multi-parameter model
		#rand_draw = random.random()
		#if rand_draw >= (1 - dup_mut) ** (max(0.5, node.dist) * gen_per_snp): #NEW: LEN 0 NODES ARE TREATED AS LEN 0.5!
		#	node.cnv_state = 'deletion'
		###
		
		###
		#Generational mutation including fractional generations
		#Works for multi-parameter model
		gens = int(node.dist * gen_per_snp)
		frac_gen = node.dist * gen_per_snp - gens #The fractional remainder of a generation
		
		for gen in xrange(gens):
			rand_draw = random.random()
			if node.cnv_state == 'reference':
				if rand_draw < dup_mut:
					node.cnv_state = 'mutation'

			else:
				if rand_draw < dup_rev:
					node.cnv_state = 'reference'
		
		rand_draw = random.random()
		if node.cnv_state == 'reference':
			if rand_draw >= (1 - dup_mut) ** frac_gen:
				node.cnv_state = 'mutation'
		else:
			if rand_draw >= (1 - dup_rev) ** frac_gen:
				node.cnv_state = 'reference'
		
		
		
		'''
		if node.cnv_state != orig_state:
			node.set_style(mutstyle)
		if node.is_leaf() and node.cnv_state != 'reference':
			node.add_face(ete3.TextFace(' %s' %(node.cnv_state[:3])), 1, 'branch-right')
			node.cnv_state = 'mutation'
		'''
	return t
	

# def haplogroup_pcts(tree):
# 	haps = {}
# 	for leaf in tree.iter_leaves():
# 		leaf_hap = leaf.name.split('_')[1][0]
# 		if leaf_hap not in haps:
# 			haps[leaf_hap] = [0,0.]
# 		haps[leaf_hap][1] += 1
# 		if leaf.cnv_state != 'reference':
# 			haps[leaf_hap][0] += 1
# 	return haps


def max_mut_age(leaf):
	if leaf.cnv_state == 'reference':
		return 0.
	age = leaf.dist
	mut_node = leaf.up
	while mut_node.cnv_state != set(['reference']):
		age += mut_node.dist
		mut_node = mut_node.up
		if mut_node.is_root():
			break
	return age

def min_mut_age(leaf):
	if leaf.cnv_state == 'reference':
		return 0.
	age = leaf.dist
	mut_node = leaf.up
	while mut_node.cnv_state != 'reference':
		age += mut_node.dist
		mut_node = mut_node.up
	return 'THIS FUNCTION NEEDS WORK'
	return age

def mut_neighbors(tree):
	leaves = tree.get_leaves()
	prev_state = 'ref' if leaves[0].cnv_state == 'reference' else 'mut'
	change_count = 0
	for leaf in leaves[1:]:
		leaf_state = 'ref' if leaf.cnv_state == 'reference' else 'mut'
		if leaf_state != prev_state:
			change_count += 1
		prev_state = leaf_state
	return change_count


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
	import argparse
	import pickle
# 	import matplotlib
# 	matplotlib.use('Agg')
# 	import matplotlib.pyplot as plt
	import os
	import numpy as np
	
	parser = argparse.ArgumentParser(description = "Simulates amplicon CNV evolution of a tree")
	parser.add_argument("-t", "--tree", default = "/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Tree_branch_adjustment/adjusted_tree_all_branches.nwk", help = "Newick file of tree")
	parser.add_argument("-u", "--dup_mut", help = "Duplication mutation rate (events per generation)", type = float, default = .00043)
	parser.add_argument("-e", "--del_mut", help = "Deletion mutation rate (events per generation)", type = float, default = .000)
	parser.add_argument("-r", "--dup_rev", help = "Duplication reversion rate (events per generation)", type = float, default = .00086)
	parser.add_argument("-v", "--del_rev", help = "Deletion reversion rate (events per generation)", type = float, default = .000)
	parser.add_argument("-i-", "--iters", help = "Number of iterations the tree is mutated", type = int, default = 1000)
	args = parser.parse_args()
	
	
	tree_structure = ete3.Tree(args.tree, format=1)
	
	yvarfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000_Ys_variant_men.txt', 'r')
	yvars = {}
	for line in yvarfile:
		if not line.startswith('#'):
			data = line.split('\t')		
			yvars[data[0]] = data[1]
		
	good_leaves = []
	for leaf in tree_structure.iter_leaves():
		leafname = leaf.name.split('_')[0]
		if leafname in yvars:
			good_leaves.append(leaf.name)
	tree_structure.prune(good_leaves)


	
	#out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/'
	#out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/Adjusted_branch_lengths/two_param_sim_stats/'
	out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/Adjusted_branch_lengths/one_param_sim_stats/'
	
	outfile = open('%sdu%s_dur%s_10000_sims.txt' %(out_dir,
											 str(args.dup_mut).split('.')[-1],
											 str(args.dup_rev).split('.')[-1]), 'w')
	
	for thing in xrange(args.iters):
		tree = evolve_tree(tree_structure, dup_mut=args.dup_mut, del_mut=args.del_mut, dup_rev=args.dup_rev, del_rev=args.del_rev)
		sim_muts = len([x for x in tree.get_leaves() if x.cnv_state != 'reference'])
		sim_fitch = fitch_up(tree)
		#sim_neighbor = mut_neighbors(tree)
		
		#fitch_down(tree)
		#sim_ages = 0
		#for leaf in tree.iter_leaves():
		#	sim_ages += max_mut_age(leaf)
		
		outfile.write('%i\t%i\t%i\t%i\n' %(sim_muts, sim_fitch, 0, 0))#sim_neighbor, sim_ages))
	
	outfile.close()
	
	
	
	