Branch Sorting Method contains:

branch_sorting.py: Calculates a p-value for likelihood of mutation distribution in the empirical 1000 Ys tree. Does so by sorting all branches of the tree from oldest to youngest, and then compares the CDF of mutation events over time to a null model where those events are expected to occur neutrally (that is, same probability at any time throughout the evolutionary history). Using a KS test between the observed and null CDFs, we get a p-value.

Requires two inputs: the phylogenetic tree of 1000 Genomes Y chromosomes (Supplemental Data 1 in Teitz et al. 2018) and a file of the 1000 Genomes males and an annotation of their amplicon variants (provided in this directory as 1000_Ys_variant_men.txt). Running this script also prints out the number of amplicon CNV mutation events that occur in the first half of the tree and the entire tree, and the sum of branch lengths in the entire tree. The is also an option (-s) to save a list of the KS test p-value of each time the branches of the tree are shuffled, so that analyses can be run on that distribution.


usage: Y_tree_evolution_pvalue.py [-h] [-t TREE] [-o OUTPUT_DIRECTORY]
                                  [-v VARIANT_ANNOTATION] [-n] [-s]

optional arguments:
  -h, --help            show this help message and exit
  -t TREE, --tree TREE  Phylogenetic tree of Y chromosomes in newick format
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Directory to write output files to (default: current
                        directory)
  -v VARIANT_ANNOTATION, --variant_annotation VARIANT_ANNOTATION
                        File containing individuals annotated with amplicon
                        variants
  -n, --no_plots        Do not create plots of real and shuffled branches
  -s, --shuffled_pvals  Create file listing p-values of 1000 shuffles
  


branch_sorting_with_simulation.py: Simulates mutation over the Y-chromosome phylogenetic tree and performs the branch sorting analysis on each simulation.


1000_Ys_variant_men.txt: File of the 1000 Genomes males and an annotation of their amplicon variants used as input for branch_sorting.py and branch_sorting_with_simulation.py.