#!/usr/bin/python

'''
Takes vectors files and makes a pandas dataframe
'''

import pandas as pd
import sys

regions = ['Haplogroup',
			'Sub-Haplogroup',
			'Single-copy_depth',
			'IR3',
			'TSPY',
			'IR1',
			'P8',
			'P7',
			'P6',
			'P5',
			'IR5',
			'P4',
			'DYZ19',
			'RBMY',
			'IR2',
			'Blue',
			'Teal',
			'Green',
			'Red',
			'DAZ',
			'Gray',
			'Yellow',
			'Ctrl_reg_0',
			'Ctrl_reg_1',
			'Ctrl_reg_2',
			'Ctrl_reg_3',
			'Ctrl_reg_4',
			'Ctrl_reg_5']

def vecfile_to_vector(vfile, name, haplogroup, depthfile=None):
	with open(vfile, 'r') as vf:
		vfdata = vf.readlines()
	if depthfile:
		with open(depthfile) as ovf:
			depth = float(ovf.readline().rstrip().split()[-1])
	else:
		depth = float(vfdata[0].rstrip().split()[-1])
	vec = [name, haplogroup, haplogroup[0], depth]
	vec += [float(x) for x in vfdata[2].rstrip().split()]
	return vec





'''

regions = ['Haplogroup', 'Sub-Haplogroup', 'Single-copy_depth']
regions += 'IR3    TSPY    IR1     P8      P7      P6      P5  IR5    P4      DYZ19   RBMY    IR2     Blue    Teal    Green   Red     DAZ     Gray Yellow Ctrl_reg_0      Ctrl_reg_1      Ctrl_reg_2      Ctrl_reg_3      Ctrl_reg_4      Ctrl_reg_5'.split()
pdvectors = pd.DataFrame([x[1:] for x in vectors], index=[x[0] for x in vectors], columns=regions)
pdvectors.to_pickle('1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls.pickle')
'''