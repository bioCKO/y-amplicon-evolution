#!/usr/bin/python

import sys
import matplotlib.pyplot as plt

insam = open(sys.argv[1], 'r')
outfile = sys.argv[2]
cutoff = int(sys.argv[3])
read = 'chrY.r0'
readcount = 0
repvec = []

repeats = open('/lab/page/lsteitz/Reference_genomes/repeat_masker_good_regions.bed', 'r')
goodlocs = set([])
for line in repeats:
	data = line.split()
	if data[0] == 'chrY':
		for base in xrange(int(data[1]),int(data[2])-10):
			goodlocs.add(base+1)

for line in insam:
	if line.startswith('@'):
		continue
	if line.split()[0] != read:
		if int(read.split('.r')[1]) in goodlocs:
			repvec.append(readcount)
		else:
			repvec.append(0)
		read = line.split()[0]
		readcount = 0
        if line.split()[1] != '4' and int(line.split()[11].split(':')[-1]) >= cutoff:
                readcount += 1

insam.close()
repvec = repvec[1:]
repvec.append(readcount)

'''
figs = 57
for fig in xrange(1,figs+1):
	plt.clf()
	plt.figure(figsize=(200,4))
	plt.plot(repvec[len(repvec)/figs*fig:len(repvec)/figs*(fig+1)])
	plt.savefig(outfile + '_%i.png' %(fig))
'''

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
