#!/usr/bin/python

'''
Takes sam file created by find_unique_regions.py and makes a bedGraph file of 
regions of the Y chromosome and how many times they map in the human genome.
'''

import sys

insam = open(sys.argv[1], 'r')
outfile = sys.argv[2]
cutoff = int(sys.argv[3])

read = ''
readcount = 0
repvec = []

for line in insam:
	if line.startswith('@'):
		continue
	if line.split()[0] != read:
		repvec.append(readcount)
		read = line.split()[0]
		readcount = 0
	if line.split()[1] != '4' and int(line.split()[11].split(':')[-1]) >= cutoff:
		readcount += 1

insam.close()
repvec = repvec[1:]
repvec.append(readcount)


outname = open('%s.bedgraph' %(outfile), 'w')	
startrep = 0
startpos = 0
for pos, rep in enumerate(repvec):
	if rep != startrep:
		if startrep != 0:
			outname.write('chrY\t%i\t%i\t%i\n' %(startpos, pos, startrep))
		startrep = rep
		startpos = pos
outname.close()
