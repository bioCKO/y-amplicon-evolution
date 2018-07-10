#!/usr/bin/python

'''
Simulation of genetic drift over a range of population sizes and mutation rates 
used to simulate haplogroup A00 amplicon mutation
'''

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

outname = sys.argv[1]
outfile = open(outname, 'w')
outfile.write('m\tN\t% Mutant\n')
#Range of mutation rates
m_range = [10**-x for x in xrange(1,7)] + [2. * 10**-x for x in xrange(1,7)] + [4. * 10**-x for x in xrange(1,7)] + [8. * 10**-x for x in xrange(1,7)] + [0]
m_range.sort()


for N in [100, 1000, 5000, 10000, 50000, 100000, 1000000]: #Range of population sizes
	for m in m_range:
		percent_mut = []
		for rep in xrange(100):
			p = 0 #starting fraction with mutations
			q = 1 #starting fraction without mutations
			#r = .0001 #term for reversion; not used in published results
			gens = 10000 #Also simulated 8000 generations; results were similar


			for generation in xrange(gens):
				p = np.random.binomial(N, p + m*q)/float(N)
				#p = np.random.binomial(N, (1-r)*p + m*q)/float(N) #Model with reversion; not used in published results
				if p == 1:
					break
				
				q = 1-p

			percent_mut.append(p)
		outfile.write('%f\t%i\t%f\n' %(m, N, np.mean(percent_mut)))
		print '%f\t%i\t%f' %(m, N, np.mean(percent_mut))
outfile.close()


