#!/usr/bin/python

'''
Get probability distribution from fragment size histogram text file
'''
import sys
import numpy as np
import random
infile = open(sys.argv[1], 'r')
hist = []

def build_hist_dist(infile):
	start = 0
	for line in infile:
		if line.startswith('#'):
			continue
		data = [int(x) for x in line.rstrip().split()]
		hist.append([data[0],data[1]])

	binsize = hist[1][0] - hist[0][0]

	windows = []
	probs = []
	for win, prob in hist:
		for size in xrange(binsize):
			windows.append(win+size)
			probs.append(prob)

	tot = float(sum(probs))
	probs = np.cumsum([x/tot for x in probs])
	
	return windows, probs
	
windows, probs = build_hist_dist(infile)
fragfake = []
for thing in xrange(24000):
	randfragsize = windows[np.searchsorted(probs, random.random())]
	fragfake.append(randfragsize)
import matplotlib.pyplot as plt
plt.clf()
plt.hist(fragfake, bins=range(min(fragfake)-10, max(fragfake)+10,max(1,(max(fragfake)-min(fragfake))/200)))
plt.savefig('FULLTEST_SIM.png')
