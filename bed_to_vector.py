#!/usr/bin/python

'''
Reads a bedGraph file of sequencing depth and returns 
a linear array of depth at each base in a chromosome.
'''

def bed_to_depthvec(bedfile, chrom, regfile=None):
	depthvec = []
	last = -1
	with open(bedfile, 'r') as infile:
		for line in infile:
			data = line.split()
			if data[0] == chrom:
				depthvec += [0] * (int(data[1]) - last - 1)
				for thing in xrange(int(data[1]), int(data[2])):
					depthvec.append(float(data[3]))
					last = thing
	if regfile:
		with open(regfile, 'r') as infile:
			goodregs = []
			for line in infile:
				data = line.split()
				goodregs.append([int(data[1]), int(data[2])])
		tempvec = []
		for reg in goodregs:
			tempvec += depthvec[reg[0]:reg[1]]
		depthvec = tempvec
	return depthvec

