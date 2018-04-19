#!/usr/bin/python

import pickle

repeats = open('/lab/page/lsteitz/Reference_genomes/repeat_masker_good_regions.bed', 'r')
goodlocs = set([])
for line in repeats:
	data = line.split()
	if data[0] == 'chrY':
		for base in xrange(int(data[1])+10,int(data[2])-10):
			goodlocs.add(base+1)

repeats.close()

outfile = open('/lab/solexa_page/lsteitz/chrY_repeatmasker_set.pickle', 'w')
pickle.dump(goodlocs, outfile)
outfile.close()
