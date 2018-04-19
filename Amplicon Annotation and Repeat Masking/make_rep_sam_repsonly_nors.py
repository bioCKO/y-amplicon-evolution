#!/usr/bin/python

import sys
import matplotlib.pyplot as plt

insam = open(sys.argv[1], 'r')
outname = sys.argv[2]
cutoff = int(sys.argv[3])

repeats = open('/lab/page/lsteitz/Reference_genomes/repeat_masker_good_regions.bed', 'r')
goodlocs = set([])
for line in repeats:
	data = line.split()
	if data[0] == 'chrY':
		for base in xrange(int(data[1]),int(data[2])-10):
			goodlocs.add(base+1)

outfile = open(outname, 'w')

for line in insam:
	if line.startswith('@'):
		outfile.write(line)
		continue
	if (line.split()[0].split('.r')[1] != line.split()[3] or line.split()[2] != 'chrY') and int(line.split()[0].split('.r')[1]) in goodlocs and line.split()[1] != '4' and int(line.split()[11].split(':')[-1]) >= cutoff:
		outfile.write(line)

outfile.close()
insam.close()
