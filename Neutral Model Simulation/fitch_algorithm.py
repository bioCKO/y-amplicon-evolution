#!/nfs/apps/anaconda2/bin/python

'''Local hashbang: #!/usr/local/bin/miniconda2/bin/python'''

'''
Run Fitch's algorithm on 1000 Ys tree with Y CNV's annotated 
and return minimum number of mutations necessary in tree
'''

import ete3

def fitch_count_events(tree):
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

if __name__ == "__main__":
	import ast
	
	#newick = '/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.nwk'
	newick = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/tree.all.nwk'
	
	#Make dict of male IDs to Y CNV states
	yvarfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000_Ys_individual_sv_file.txt', 'r')
	yvars = {}
	ref_cns = (2, 2, 2, 2, 2, 4, 2, 2, 4, 2, 3, 2, 2) #Reference copy numbers for use when counting only dels or dups
	for line in yvarfile:
		if not line.startswith('#'):
			data = line.split('\t')		
			yvars[data[0]] = data[1]
			#THESE TWO LINES MAKES IT ONLY COUNT AZFc MUTATIONS:
			#if not data[1] == 'v1' and data[2].rstrip().endswith('4, 2, 3, 2, 2)'):
			#	yvars[data[0]] = 'v1'
			#THESE TWO LINES MAKES IT ONLY COUNT DELETIONS (>=) OR DUPLICATIONS (<=):
			#if any([ast.literal_eval(data[2].rstrip())[x] > ref_cns[x] for x in xrange(len(ref_cns))]):
			#	yvars[data[0]] = 'v1'


	t = ete3.Tree(newick, format=9)
	
	#Annotate each leaf's CNV state and remove nodes not in my data
	good_leaves = []
	for leaf in t.iter_leaves():
		leafname = leaf.name.split('_')[0]
		if leafname in yvars:
			good_leaves.append(leaf.name)
			leaf.add_feature('cnv_state', set([yvars[leafname]]))
	t.prune(good_leaves)
	
	
	#Traverse tree and run first part of Fitch's algorithm

	event_count = fitch_count_events(t)

	print len([x for x in t.get_leaves() if x.cnv_state != set(['v1'])])
	print event_count