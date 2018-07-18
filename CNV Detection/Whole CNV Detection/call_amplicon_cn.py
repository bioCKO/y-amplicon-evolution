#!/usr/bin/python

'''
Calls amplicon copy number from normalized amplicon depth
'''

import pickle
import argparse
import numpy as np
import pandas as pd

def call_cn(ampname, ampdepth):

	cn_2 = ['IR3', 'IR1', 'P8', 'P7', 'P6', 'P5', 'P4', 'IR2', 'Teal', 'Gray', 'Yellow']
	cn_3 = ['Green']
	cn_4 = ['IR5', 'Blue', 'Red']

	copy_numbers = {}
	for amp in cn_2:
		copy_numbers[amp] = 2
	for amp in cn_3:
		copy_numbers[amp] = 3
	for amp in cn_4:
		copy_numbers[amp] = 4
			
	if np.isnan(ampdepth) or np.isinf(ampdepth):
		return 0
	else:
		return int(ampdepth * copy_numbers[ampname] + 0.5)
		

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Calls amplicon copy number from normalized amplicon depth")
	
	parser.add_argument("infile", help = "pickled pandas dataframe of amplicon depths")
	parser.add_argument("-o", "--outfile", default = "./amplicon_calls.pickle", help = "Output filename")
	
	args = parser.parse_args()

	complex_and_control = ['TSPY', 'DYZ19', 'RBMY', 'DAZ', 
						   'Ctrl_reg_0', 'Ctrl_reg_1', 
						   'Ctrl_reg_2', 'Ctrl_reg_3', 'Ctrl_reg_4', 
						   'Ctrl_reg_5']

	vectors = pd.read_pickle(args.infile)
	vectors = vectors.filter([x for x in vectors.columns if x not in complex_and_control])

	for ind in vectors.index:
		for amp in vectors.columns[3:]:
			norm_cov = vectors.loc[ind, amp]
			vectors.set_value(ind, amp, call_cn(amp, norm_cov))

	outfile = open('%s' %(args.outfile), 'w')
	pickle.dump(vectors, outfile)
	outfile.close()

