#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

percent_mut = []

for rep in xrange(1000):
	p = 0
	q = 1
	m = .00043
	r = .0001
	N = 10000
	gens = 10000

	results = []

	for generation in xrange(gens):
		p = np.random.binomial(N, p + m*q)/float(N)
		#p = np.random.binomial(N, (1-r)*p + m*q)/float(N) #With Reversion
		q = 1-p
		results.append(p)

	plt.plot(results)
	percent_mut.append(p)
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/drift/m%s_N%i_%igens_rev.png' %(str(m)[1:], N, gens))

plt.clf()
plt.hist(percent_mut)
plt.title('Mean = %f' %(np.mean(percent_mut)))
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/drift/m%s_N%i_%igens_hist_rev.png' %(str(m)[2:], N, gens))


