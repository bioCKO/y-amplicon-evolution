#!/nfs/apps/anaconda2/bin/python

'''
Gets a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree.
Does so by sorting all branches of the tree from oldest to youngest, and then compares 
the CDF of mutation events over time to a null model where those events are expected to 
occur neutrally (that is, same probability at any time throughout the evolutionary history).
Using a KS test between the observed and null CDFs, we get a p-value.
'''

import ete3
import pickle
import random
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from scipy.stats import kstest
from Y_tree_evolution_sim_4 import fitch_up, fitch_down

out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/no_one_amps'

tree = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.branch.lengths.nwk', format=1)


#Annotate tree with variants
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

#No partial support
one_amp_list = ['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA19774', 'NA18982', 'NA19198', 'NA19451', 'HG04219', 'HG03021', 'HG03015', 'NA18977']

#All
#one_amp_list['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG02512', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA21093', 'HG01842', 'NA21112', 'HG01767', 'HG01530', 'NA21111', 'NA19774', 'NA18982', 'NA19198', 'NA12748', 'HG03660', 'NA19451', 'HG04219', 'NA20884', 'HG03021', 'HG03015', 'NA18977', 'NA20536']


for leaf in tree.iter_leaves():
	leafname = leaf.name.split('_')[0]
	if leafname in yvars:
		if yvars[leafname] == 'v1' or leafname in one_amp_list:
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


'''
no_states = []
ctcount = 0
for branch in tree.traverse():
	no_states.append([branch.cnv_state, len(branch.cnv_state)])
	if len(branch.cnv_state) > 1:
		for node in branch.traverse():
			if node.is_leaf():
				nodename = node.name.split('_')[0]
				node.add_face(ete3.TextFace(' %s' %(yvars[nodename])), 1, 'branch-right')
			else:
				node.add_face(ete3.TextFace(node.name), 1, 'branch-right')

	
		bl_style = ete3.TreeStyle()
		bl_style.show_branch_length = True
		out_dir = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/'
		ctcount += 1
		branch.render('%scomplex_trees_%i.pdf' %(out_dir, ctcount), tree_style=bl_style)
'''

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
						 #'HG01097_G2', #remove for one amp filtering
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
						 #'NA18982_O3', #remove for one amp filtering
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

tot_age = sum([x[1] for x in branches])
first_half_muts = 0
age = 0
branch_idx = 0
while age < tot_age/2:
	if branches[branch_idx][2] == 'yes':
		first_half_muts += 1
	age += branches[branch_idx][1]
	branch_idx += 1

print first_half_muts

while age < tot_age:
	if branches[branch_idx][2] == 'yes':
		first_half_muts += 1
	age += branches[branch_idx][1]
	branch_idx += 1

print first_half_muts



plt.figure(figsize=[10,4])

age = 0
for branch in branches:                                                                              
	plt.plot([age, age + branch[1]], [.2, .2])
	if branch[3] == 1:
		plt.plot([age], [0.3], 'b.')
	age += branch[1]
	if branch[2]  == 'yes':
		plt.plot([age], [0.25], 'r.')

plt.axis([0,age,0.1,0.4])

plt.savefig('%s/mutations_from_leaves_sorted_halfway.png' %(out_dir), dpi=1000)


plt.clf()

plt.figure(figsize=[7,7])
age = 0
mut_times = []
for branch in branches:
	if branch[2] == 'yes':
		mut_times.append(age + branch[1]/2.)
	age += branch[1]

ks_p = kstest(np.array(mut_times)/float(age), 'uniform')
print ks_p
print age
print fitch_events
print fitch_events/float(age)


plt.plot(mut_times, np.arange(0,1,1./len(mut_times)))
plt.plot([0,age],[0,1], 'r')
plt.axis([0, age, 0, 1])
plt.xlabel('Evolutionary time\n(SNPs)')
plt.ylabel('Fraction of variants observed')
plt.title(r'$p = 1.140\times10^{-9}$ (KS test)')
plt.savefig('%s/mutations_over_time.png' %(out_dir), dpi=300)

	
#Make figure of shuffled branches
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
	
	


	if subfig in [1, 10, 100, 500]:
		continue
		plt.plot([0,age],[0,1], 'r', lw=2)
		plt.axis([0, age, 0, 1])
		plt.savefig('%s/mutations_over_time_shuffled_%ishufs.png' %(out_dir, subfig), dpi=300)
	
	
	
	
plt.plot([0,age],[0,1], 'r', lw=2)
plt.axis([0, age, 0, 1])
#plt.savefig('%s/mutations_over_time_shuffled.png' %(out_dir), dpi=300)


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
plt.xlabel('Evolutionary time\n(SNPs)')
plt.ylabel('Fraction of variants observed')
#plt.title(r'$p = 1.140\times10^{-9}$ (KS test)')
plt.title(r'$p = %s\times10^{-%i}$ (KS test)' %(('%3e' %ks_p[1])[:5], int(str(ks_p[1]).split('-')[-1])))

plt.savefig('/%s/mutations_over_time_with_shuffled_high_confidence_calls.png' %(out_dir), dpi=300)
plt.savefig('/%s/mutations_over_time_with_shuffled_high_confidence_calls.svg' %(out_dir))


# import pickle
# outfile = open('%s/shuffled_ks_pvals.pickle' %(out_dir), 'w')
# pickle.dump(ks_ps, outfile)
# outfile.close()
'''
#Make combined figure
plt.clf()

plt.figure(figsize=[7,7])

for subfig in xrange(10000):
	random.shuffle(branches)
	age = 0
	shuf_mut_times = []
	for branch in branches:
		if branch[2] == 'yes':
			shuf_mut_times.append(age + branch[1]/2.)
		age += branch[1]

	plt.plot(shuf_mut_times, np.arange(0,1,1./len(shuf_mut_times)), 'gray', alpha=0.3)

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
plt.xlabel('Evolutionary time\n(SNPs)')
plt.ylabel('Fraction of variants observed')
plt.title(r'$p = 3.466\times10^{-10}$ (KS test)')
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Y_tree_pval/mutations_over_time_with_sims.png')
'''	






