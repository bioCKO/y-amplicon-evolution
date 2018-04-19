#!/nfs/apps/anaconda2/bin/python


import ete3


rstyle = ete3.NodeStyle()
rstyle['size'] = 0
allstyle = ete3.TreeStyle()
allstyle.show_leaf_name = False

tree_s = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/a1.names.nwk', format=1)
rescue = tree_s&"NA18960_D2"
rescue = rescue.up

rescued_man = 'NA18960'

for node in rescue.traverse():
    node.set_style(rstyle)
for node in rescue.iter_leaves():
    node.name = node.name.split('_')[0]
for node in rescue.iter_leaves():
    if node.name == rescued_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='green'), column=0)
    else:
        node.add_face(ete3.TextFace(node.name, fgcolor='red'), column=0)

rescue.render('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Rescued deletions/D2.svg', tree_style=allstyle)

tree_s = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/b.names.nwk', format=1)
rescue = tree_s&"HG03646_H1"
rescue = rescue.up
var_man = 'HG03646'
rescued_man = 'HG03745'

for node in rescue.traverse():
    node.set_style(rstyle)
for node in rescue.iter_leaves():
    node.name = node.name.split('_')[0]
for node in rescue.iter_leaves():
    if node.name == rescued_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='green'), column=0)
    elif node.name == var_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='red'), column=0)
    else:
        node.add_face(ete3.TextFace(node.name, fgcolor='black'), column=0)
rescue.render('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Rescued deletions/H1.svg', tree_style=allstyle)

tree_s = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/e1.names.nwk', format=1)
rescue = tree_s&"HG02696_Q1a"
rescue = rescue.up
rescued_man = ['HG02696', 'HG02134']
var_man = ['HG02116', 'HG01944']
for node in rescue.traverse():
    node.set_style(rstyle)
for node in rescue.iter_leaves():
    node.name = node.name.split('_')[0]
for node in rescue.iter_leaves():
    if node.name in rescued_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='green'), column=0)
    elif node.name in var_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='red'), column=0)
rescue.render('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Rescued deletions/Q1a.svg', tree_style=allstyle)

tree_s = ete3.Tree('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/1000Y.trees/snp.branch.lengths.fixed/d1.names.nwk', format=1)
rescue = tree_s&"HG00351_N1"
rescue = rescue.up
rescue = rescue.up
rescued_man = 'HG00280'
for node in rescue.traverse():
    node.set_style(rstyle)
for node in rescue.iter_leaves():
    node.name = node.name.split('_')[0]
for node in rescue.iter_leaves():
    if node.name == rescued_man:
        node.add_face(ete3.TextFace(node.name, fgcolor='green'), column=0)
    else:                  
        node.add_face(ete3.TextFace(node.name, fgcolor='red'), column=0)
rescue.render('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Rescued deletions/N1.svg', tree_style=allstyle)
