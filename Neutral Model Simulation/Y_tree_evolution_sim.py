#!/nfs/apps/anaconda2/bin/python

'''
Hashbang for running on local machine: #!/usr/local/bin/miniconda2/bin/python
'''
'''
Simulates amplicon mutation in 1000 Ys tree to get haplogroup CNV % distributions
'''

import ete3
import random
from fitch_algorithm import fitch_count_events
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
		for generation in xrange(int(node.dist * gen_per_snp)): #maybe change this to a power later--would require editing evolve_node function
			evolve_node(node, dup_mut=dup_mut, del_mut=del_mut, dup_rev=dup_rev, del_rev=del_rev)
		if node.cnv_state != orig_state:
			node.set_style(mutstyle)
		if node.is_leaf() and node.cnv_state != 'reference':
			node.add_face(ete3.TextFace(' %s' %(node.cnv_state[:3])), 1, 'branch-right')
	
	return t
	

def haplogroup_pcts(tree):
	haps = {}
	for leaf in tree.iter_leaves():
		leaf_hap = leaf.name.split('_')[1][0]
		if leaf_hap not in haps:
			haps[leaf_hap] = [0,0.]
		haps[leaf_hap][1] += 1
		if leaf.cnv_state != 'reference':
			haps[leaf_hap][0] += 1
	return haps


if __name__ == '__main__':
	import argparse
	import pickle
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	import os
	import numpy as np
	
	parser = argparse.ArgumentParser(description = "Simulates amplicon CNV evolution of a tree")
	parser.add_argument("-t", "--tree", default = "/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.branch.lengths.nwk", help = "Newick file of tree")
	parser.add_argument("-u", "--dup_mut", help = "Duplication mutation rate (events per generation)", type = float, default = .00043)
	parser.add_argument("-e", "--del_mut", help = "Deletion mutation rate (events per generation)", type = float, default = .000)
	parser.add_argument("-r", "--dup_rev", help = "Duplication reversion rate (events per generation)", type = float, default = .00086)
	parser.add_argument("-v", "--del_rev", help = "Deletion reversion rate (events per generation)", type = float, default = .000)
	parser.add_argument("-i-", "--iters", help = "Number of iterations the tree is mutated", type = int, default = 1000)
	args = parser.parse_args()
	
	# out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/all'
# 	#out_dir = '/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/all'
# 	
# 	if args.dup_mut > 0:
# 		out_dir += '_du' + str(args.dup_mut).split('.')[-1]
# 	if args.dup_rev > 0:
# 		out_dir += '_dur' + str(args.dup_rev).split('.')[-1]
# 	if args.del_mut > 0:
# 		out_dir += '_de' + str(args.del_mut).split('.')[-1]
# 	if args.del_rev > 0:
# 		out_dir += '_der' + str(args.del_rev).split('.')[-1]
# 	out_dir += '/'
# 		
# 	if not os.path.exists(out_dir):
# 		os.mkdir(out_dir)
	
	tree_structure = ete3.Tree(args.tree, format=1)
	
	ts = ete3.TreeStyle()
	ts.show_branch_length = True
	
	#haplogroups = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'L', 'N', 'O', 'Q', 'R', 'T']
	#Above: all haplogroups; below: haplogroups with >10 men
	haplogroups = ['C', 'D', 'E', 'G', 'H', 'I', 'J', 'L', 'N', 'O', 'Q', 'R']
	happcts = {x:[] for x in haplogroups}
	happctvecs = []
	
	for thing in xrange(args.iters):
		tree = evolve_tree(tree_structure, dup_mut=args.dup_mut, del_mut=args.del_mut, dup_rev=args.dup_rev, del_rev=args.del_rev)
		# print [x.cnv_state for x in tree.iter_leaves()].count('duplication')/float(len(tree))
# 		print fitch_count_events(tree)
# 		print  len([x for x in tree.get_leaves() if x.cnv_state != 'reference'])
		tree_happcts = haplogroup_pcts(tree)
		tree_happcts = {x : tree_happcts[x] for x in haplogroups}
		tree_fitch_events = fitch_count_events(tree)
		happctvecs.append([tree_happcts, tree_fitch_events])
		for hap in happcts:
			happcts[hap].append(tree_happcts[hap][0]/tree_happcts[hap][1])
		#tree.render('/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/all_du00043_dur00086/all_du00043_dur00086_%i.pdf' %(thing), tree_style=ts)
		# if thing%25 == 0:
		#	print thing
	# outhappcts = open('/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/all_du00043_dur00086/happcts_2.list', 'w')
# 	pickle.dump(happcts, outhappcts)
# 	outhappcts.close()
# 	outhappctsvecs = open('/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/all_du00043_dur00086/happctvecs.list', 'w')
# 	pickle.dump(happctvecs, outhappctsvecs)
# 	outhappctsvecs.close()
	
	#realpcts = {'A': 0.1428571428571429, 'C': 0.1785714285714286, 'B': 0.0, 'E': 0.1149068322981367, 'D': 1.0, 'G': 0.2777777777777778, 'F': 0.0, 'I': 0.0625, 'H': 0.16129032258064513, 'J': 0.1578947368421053, 'L': 0.2592592592592593, 'O': 0.20812182741116747, 'N': 0.967741935483871, 'Q': 0.11111111111111116, 'R': 0.06687898089171973, 'T': 0.375}
	#Above: all haplogroups; below: haplogroups with >10 men
	realpcts = {'C': 0.1785714285714286,'E': 0.1149068322981367, 'D': 1.0, 'G': 0.2777777777777778, 'I': 0.0625, 'H': 0.16129032258064513, 'J': 0.1578947368421053, 'L': 0.2592592592592593, 'O': 0.20812182741116747, 'N': 0.967741935483871, 'Q': 0.11111111111111116, 'R': 0.06687898089171973}
	
	ordered_hpv = [[] for x in haplogroups]
	for thing in happctvecs:
		hpv_vec = sorted([x[0]/x[1] for x in thing[0].values()], reverse=True)
		for window in xrange(len(ordered_hpv)):
			ordered_hpv[window].append(hpv_vec[window])
	
	
	#Test to see how like or unlike the real data is to the simulated data.
	#For each haplogroup category (highest variant pct, 2nd highest, etc.),
	#take the real value and see how what fraction of the simulated values in 
	#that category fall within 0.1 of the real value. Multiply each of those 
	#together to get a single joint value. Then compare that joint value to the 
	#joint values of each simulated run and see where it falls in that distribution.
	real_joint_value = 1
	for pct_idx, pct in enumerate(sorted(realpcts.values(), reverse=True)):
		#How many sims are within 10% of data
		#within_tenpct = len(filter(lambda x: abs(x - pct) <= 0.1, ordered_hpv[pct_idx]))
		#real_joint_value *= (within_tenpct + 1)/(float(args.iters) + 1) #The +1s are to prevent a single 0 from wiping out the entire joint product
		#How many sims are between data and median of sims
		hpv_median = np.median(ordered_hpv[pct_idx])
		if pct < hpv_median:
			off_med = len(filter(lambda x: pct < x < hpv_median, ordered_hpv[pct_idx]))
		else:
			off_med = len(filter(lambda x: pct > x > hpv_median, ordered_hpv[pct_idx]))
		real_joint_value *= (args.iters/2 - off_med + 1)/(float(args.iters)/2 + 1)
		
	sim_joint_values = []
	
	less_likely_fitches = []
	
	for simmedpct in happctvecs:
		sim_joint_value = 1
		simmedpct_vals = sorted([x[0]/x[1] for x in simmedpct[0].values()], reverse=True)
		for pct_idx, pct in enumerate(simmedpct_vals):
			#How many sims are within 10% of data
			#within_tenpct = len(filter(lambda x: abs(x - pct) <= 0.1, ordered_hpv[pct_idx]))
			#sim_joint_value *= (within_tenpct + 1)/(float(args.iters) + 1) #The +1s are to prevent a single 0 from wiping out the entire joint product
			hpv_median = np.median(ordered_hpv[pct_idx])
			if pct < hpv_median:
				off_med = len(filter(lambda x: pct < x < hpv_median, ordered_hpv[pct_idx]))
			else:
				off_med = len(filter(lambda x: pct > x > hpv_median, ordered_hpv[pct_idx]))
			sim_joint_value *= (args.iters/2 - off_med + 1)/(float(args.iters)/2 + 1)
		sim_joint_values.append(sim_joint_value)
		if sim_joint_value < real_joint_value:
			less_likely_fitches.append(simmedpct[1])
	
	print real_joint_value
	p_value = len(filter(lambda x: x < real_joint_value, sim_joint_values))/float(args.iters)
	print p_value
	print np.mean(less_likely_fitches), np.std(less_likely_fitches), max(less_likely_fitches), min(less_likely_fitches)
	print 'Between median and value'#'Range: 0.1' #Line to describe method and parameter(s) of probability calculation
	quit()
	
	plt.clf()
	plt.hist(sim_joint_values, bins=np.arange(min(sim_joint_values), max(sim_joint_values), (max(sim_joint_values) - min(sim_joint_values))/100))
	plt.savefig('%ssimmed_pval_dist.png' %(out_dir))
	quit()
	
	out_name = 'sim_vs_real_violin'
	if os.path.isfile('%s%s.png' %(out_dir, out_name)):
		run_label = 2
		while os.path.isfile('%s%s_%i.png' %(out_dir, out_name, run_label)):
			run_label += 1
		out_name = '%s_%i.png' %(out_name, run_label)
	
	plt.clf()
	plt.violinplot([happcts[x] for x in sorted(happcts.keys())])
	plt.plot(range(1,len(happcts)+1), [realpcts[x] for x in sorted(happcts.keys())], 'ro')
	plt.xticks(range(1,len(happcts)+1),  sorted(happcts.keys()))
	plt.savefig('%s%s.png' %(out_dir, out_name))
	
	plt.clf()
	plt.violinplot(ordered_hpv)
	plt.plot(range(1,len(haplogroups)+1), sorted(realpcts.values(), reverse=True), 'ro')
	plt.savefig('%s%s_ordered.png' %(out_dir, out_name))
	plt.clf()
	
	
	
	
	
	
		
		

	







