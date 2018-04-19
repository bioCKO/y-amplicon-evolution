#!/usr/bin/python

from Bio import SeqIO
import sys

win_size = int(sys.argv[1])

chrY_it = SeqIO.parse('/lab/page/lsteitz/Simulations/Reference_chromosomes/chrY.fa', 'fasta')
chrY = chrY_it.next().seq
read_out = 'all_chrY_reads_%ibp.fa' %(win_size)

with open(read_out, 'w') as outfile:
	readcount = 1
	for read in xrange(len(chrY) - win_size + 1):
		outfile.write('>chrY.r%i\n' %(readcount))
		outfile.write(str(chrY[read:read + int(win_size)]) + '\n')
		readcount += 1

