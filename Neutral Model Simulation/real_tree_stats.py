#!/nfs/apps/anaconda2/bin/python

'''
Calculate total mutants, fitch events, Helen's neighbor score, 
and total mutant age in real 1000 Genomes tree.
'''

import ete3
from Y_tree_evolution_sim_3 import max_mut_age, mut_neighbors, fitch_up, fitch_down
import numpy as np

tree = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.branch.lengths.nwk', format=1)

yvarfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000_Ys_variant_men.txt', 'r')
yvars = {}
for line in yvarfile:
	if not line.startswith('#'):
		data = line.split('\t')		
		yvars[data[0]] = data[1]

ref_cn = np.array((2, 2, 2, 2, 2, 4, 2, 2, 4, 2, 3, 2, 2))
dels = ['v2', 'v8', 'v11', 'v19', 
		'v21', 'v22', 'v23', 'v24', 
		'v25', 'v26', 'v37', 'v40', 
		'v41', 'v42', 'v46', 'v50', 
		'v52', 'v54']

#No partial support
one_amp_list = ['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA19774', 'NA18982', 'NA19198', 'NA19451', 'HG04219', 'HG03021', 'HG03015', 'NA18977']

#All
#one_amp_list['HG02020', 'HG02401', 'HG01031', 'NA18953', 'HG01377', 'HG00982', 'HG02512', 'HG00707', 'HG04131', 'HG01097', 'HG02390', 'HG00742', 'HG03767', 'HG02774', 'HG01028', 'NA21093', 'HG01842', 'NA21112', 'HG01767', 'HG01530', 'NA21111', 'NA19774', 'NA18982', 'NA19198', 'NA12748', 'HG03660', 'NA19451', 'HG04219', 'NA20884', 'HG03021', 'HG03015', 'NA18977', 'NA20536']

	
good_leaves = []
for leaf in tree.iter_leaves():
	leafname = leaf.name.split('_')[0]
	if leafname in yvars:
		good_leaves.append(leaf.name)
		if yvars[leafname] == 'v1' or leafname in one_amp_list:
		#if yvars[leafname] == 'v1' or yvars[leafname] in dels:
			leaf.add_feature('cnv_state', 'reference')
		else:
			#leaf.add_feature('cnv_state', yvars[leafname])
			leaf.add_feature('cnv_state', 'deldup')
tree.prune(good_leaves)



mutants = len([x for x in tree.get_leaves() if x.cnv_state != 'reference'])
fitch_events = fitch_up(tree)
neighbor_score = mut_neighbors(tree)

fitch_down(tree)
total_mut_age = 0
for leaf in tree.iter_leaves():
	total_mut_age += max_mut_age(leaf)

print 'Mutants: %i' %(mutants)
print 'Fitch events: %i' %(fitch_events)
print 'Neighbor score: %i' %(neighbor_score)
print 'Mutants\' age: %i' %(total_mut_age)
