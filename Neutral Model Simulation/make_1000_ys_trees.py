#!/usr/local/bin/miniconda2/bin/python

import ete3
import glob

#newicks = glob.glob('/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/*_EDITED.names.nwk')
tree_dir = '/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/'
newicks = ['a1.names.nwk', 'a2_EDITED.names.nwk', 'b_EDITED.names.nwk', 
		   'c.names.nwk', 'd1.names.nwk', 'd2_EDITED.names.nwk', 
		   'e1.names.nwk', 'e2.names.nwk']

#newick_dir = '/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/'
#newicks = ['tree.all.nwk']

yvarfile = open('/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000_Ys_variant_men.txt', 'r')
yvars = {}
for line in yvarfile:
	if not line.startswith('#'):
		data = line.split()
		yvars[data[0]] = data[1]

varstyle = ete3.NodeStyle()
varstyle['fgcolor'] = 'red'

for newick in newicks:
	t = ete3.Tree('%s%s' %(tree_dir, newick), format=1)
	
	#Remove nodes not in my data
	good_leaves = []
	for leaf in t.iter_leaves():
		if leaf.name.split('_')[0] in yvars:
			good_leaves.append(leaf.name)
	t.prune(good_leaves)

# 	for node in t.search_nodes():
# 		nodename = node.name.split('_')[0]
# 		if node.is_leaf():
# 			if yvars[nodename] != 'v1':
# 				node.set_style(varstyle)
# 				node.add_face(ete3.TextFace(' %s' %(yvars[nodename])), 1, 'branch-right')

	for node in t.iter_leaves():
		nodename = node.name.split('_')[0]
		if yvars[nodename] != 'v1':
			node.add_face(ete3.TextFace(' %s' %(yvars[nodename])), 1, 'branch-right')
			node.add_feature('cnv_state', True)
		else:
			node.add_feature('cnv_state', False)
				
	
	for node in t.get_monophyletic(values=[True], target_attr='cnv_state'):
		for subnode in node.traverse():
			subnode.set_style(varstyle)
					
# 	for node in t.search_nodes():
# 		nodename = node.name.split('_')[0]
# 		if node.is_leaf():
# 			if not nodename in yvars:
# 				node.detach()
# 			elif yvars[nodename] != 'v1':
# 				node.set_style(varstyle)
# 				node.add_face(ete3.TextFace(' %s' %(yvars[nodename])), 1, 'branch-right')
# 	
# 	for iteration in xrange(20):
# 		for node in t.traverse():
# 			colors = []
# 			for subnode in node.iter_descendants():
# 				colors.append(subnode.img_style['fgcolor'])
# 			if len(colors) > 0 and all([y == 'red' for y in colors]):
# 				node.set_style(varstyle)
# 	
# 	for node in t.traverse('postorder'):
# 		colors = []
# 		for subnode in node.iter_descendants():
# 			colors.append(subnode.img_style['fgcolor'])
# 		if len(colors) > 0 and all([y == 'red' for y in colors]):
# 			node.set_style(varstyle)
				
	#print len(t.get_leaves())
	
	bl_style = ete3.TreeStyle()
	bl_style.show_branch_length = True
	#t.render('%s.png' %(newick))
	out_dir = '/Volumes/page_lab/users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Final_Trees'
	t.render('%s/%s.pdf' %(out_dir, newick.split('/')[-1].split('.')[0].split('_')[0]), tree_style=bl_style)
	