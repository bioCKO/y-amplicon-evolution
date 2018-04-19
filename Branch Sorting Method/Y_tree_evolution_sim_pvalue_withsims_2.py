#!/nfs/apps/anaconda2/bin/python

'''
Gets a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree.
Does so by sorting all branches of the tree from oldest to youngest, and then compares 
the CDF of mutation events over time to a null model where those events are expected to 
occur neutrally (that is, same probability at any time throughout the evolutionary history).
Using a KS test between the observed and null CDFs, we get a p-value.
'''

import sys
import ete3
import copy
import pickle
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from scipy.stats import kstest
from Y_tree_evolution_sim_4 import fitch_up, fitch_down

def evolve_node(node, mut_rate=0.00045):
	rand_draw = random.random()
	if rand_draw < mut_rate:
		if node.mutation == 'no':
			node.mutation = 'yes'
		else:
			node.mutation = node.mutation + 'yes'

def evolve_tree(newick, mut_rate=0.00045):
	t = copy.deepcopy(newick)

	gen_per_snp = 1 / (25 * .0076) #Generations per unit of branch length/SNP; .0076 is SNPs per year, 25 is years per generation
	for node in t.traverse():
		node.add_feature('mutation', 'no')
		'''
		#Full generational simulation; use when multiple mutations per branch are allowed
		for generation in xrange(int(node.dist * gen_per_snp)): #maybe change this to a power later--would require editing evolve_node function
			evolve_node(node, mut_rate=mut_rate)
		'''
		#Shortcut simulation; more accurate times since generations aren't rounded down
		rand_draw = random.random()
		if rand_draw >= (1 - mut_rate) ** (node.dist * gen_per_snp):
			node.mutation = 'yes'
		
	return t


tree = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.branch.lengths.nwk', format=1)

yvarfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000_Ys_variant_men.txt', 'r')
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


mut_rate = float(sys.argv[1])

for sim in xrange(1000):
	sim_tree = evolve_tree(tree, mut_rate=mut_rate)
	
	branches = []
	for branch in sim_tree.traverse():
	
		distance = 0
		nextdist = branch.dist
		mutation = branch.mutation
	
		if branch.is_leaf():
			leaf_status = 1
		else:
			leaf_status = 0
	
	#	# distance from root
	# 	while not branch.is_root():
	# 		distance += branch.dist
	# 		branch = branch.up
	
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
	
	#Only a single mutation can occur per branch
	mut_events.append(len([x[2] for x in branches if x[2].startswith('yes')]))
	
	'''
	#Multiple mutations can occur per branch
	mut_ecount = 0
	for x in branches:
		if x[2].startswith('yes'):
			mut_ecount += len(x[2])/3
	mut_events.append(mut_ecount)
	'''
	
	age = 0
	mut_times = []
	for branch in branches:
		#Only a single mutation can occur per branch
		if branch[2] == 'yes':
			mut_times.append(age + branch[1]/2.)
		
		'''
		#Multiple mutations can occur per branch
		if branch[2].startswith('yes'):
			for mut_e in xrange(len(branch[2])/3):
				mut_times.append(age + branch[1]/2.)
		'''
		
		age += branch[1]
		
	
	if len(mut_times) > 0:
		ks_ps.append(kstest(np.array(mut_times)/float(age), 'uniform'))
	else:
		ks_ps.append((0,1))
	

	if len(mut_times) > 0:
		plt.plot(mut_times, np.arange(1, len(mut_times) + 1)/float(len(mut_times)), 'gray', alpha=0.3)

	if sim in [0, 9, 99, 499]:
		plt.plot([0,age],[0,1], 'r', lw=2)
		plt.axis([0, age, 0, 1])
		#plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/simmed_%isims.png' %(sim + 1))
		


plt.plot([0,age],[0,1], 'r', lw=2)
plt.axis([0, age, 0, 1])

#plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/simmed.png')



for leaf in tree.iter_leaves():
	leafname = leaf.name.split('_')[0]
	if leafname in yvars:
		if yvars[leafname] == 'v1':
			leaf.add_feature('cnv_state', set(['reference']))
		else:
			leaf.add_feature('cnv_state', set([yvars[leafname]]))

#Name unnamed nodes for ease of annotation later
branch_no = 1
for branch in tree.traverse():
	if branch.name == '':
		branch.name = str(branch_no)
		branch_no += 1

fitch_events = fitch_up(tree)
fitch_down(tree)

#Annotate branches in which mutation occurred:
for branch in tree.traverse():
	#Unambiguous mutations
	if branch.is_root():
		branch.add_feature('mutation', 'no')
		continue
	
	if len(branch.cnv_state) == 1 and len((branch.up).cnv_state) == 1 and branch.cnv_state != (branch.up).cnv_state:
		branch.add_feature('mutation', 'yes')
	
	#This is the list of manually curated complex parts of the tree, 
	#based on my knowledge of the actual variant structures.
	elif branch.name in ['HG03812_H', 
						 '20|94',
						 '266|N-M231|175',
						 '257|77',
						 '246|26',
						 '246|16',
						 'HG02696_Q1a',
						 '143|8',
						 'HG01097_G2',
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
						 'NA18982_O3',
						 'HG02410_O2',
						 'HG02408_O2',
						 'HG02394_O2',
						 'HG03196_E1b',
						 'NA18504_E1b']:
	
		branch.add_feature('mutation', 'yes')
	
	else:
		branch.add_feature('mutation', 'no')
		

branches = []
for branch in tree.traverse():
	
	distance = 0
	nextdist = branch.dist
	mutation = branch.mutation
	
	if branch.is_leaf():
		leaf_status = 1
	else:
		leaf_status = 0
	
#	# distance from root
# 	while not branch.is_root():
# 		distance += branch.dist
# 		branch = branch.up
	
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

age = 0
mut_times = []
for branch in branches:
	if branch[2] == 'yes':
		mut_times.append(age + branch[1]/2.)
	age += branch[1]

plt.plot(mut_times, np.arange(0,1,1./len(mut_times)), 'b', lw=2)
plt.plot([0,age],[0,1], 'r', lw=2)
plt.axis([0, age, 0, 1])
plt.title(r'$p = 1.140\times10^{-9}$ (KS test)')
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/simmed_with_real.svg')





#print mut_events
#print ks_ps
'''
plt.clf()
plt.figure()
plt.hist([x[1] for x in ks_ps], bins=20)
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/simmed_ks_pvals_multi_mut%f.png' %(mut_rate))

tiles = np.percentile(ks_ps, np.arange(0,100,5))
print tiles
print max(ks_ps)
print min(ks_ps)
print sorted(ks_ps)[:20]
'''

import pickle
outfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/simmed_ks_pvals.pickle', 'w')
pickle.dump(ks_ps, outfile)
outfile.close()