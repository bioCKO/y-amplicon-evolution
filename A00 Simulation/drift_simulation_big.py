#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


outfile = open('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/drift/stats_8000_gens.txt', 'w')
outfile.write('m\tN\t% Mutant\n')
m_range = [10**-x for x in xrange(1,7)] + [2. * 10**-x for x in xrange(1,7)] + [4. * 10**-x for x in xrange(1,7)] + [8. * 10**-x for x in xrange(1,7)] + [0]
m_range.sort()


for N in [100, 1000, 5000, 10000, 50000, 100000, 1000000]:
	for m in m_range:
		percent_mut = []
		for rep in xrange(100):
			p = 0
			q = 1
			#r = .0001
			gens = 8000


			for generation in xrange(gens):
				p = np.random.binomial(N, p + m*q)/float(N)
				if p == 1:
					break
				#p = np.random.binomial(N, (1-r)*p + m*q)/float(N) #With Reversion
				q = 1-p

			percent_mut.append(p)
		outfile.write('%f\t%i\t%f\n' %(m, N, np.mean(percent_mut)))
		print '%f\t%i\t%f' %(m, N, np.mean(percent_mut))
outfile.close()


