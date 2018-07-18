#!/usr/bin/python

'''
Takes individual amplicon depth files and combines then into a pandas dataframe
'''

import pandas as pd
import argparse

def depthfile_to_dfrow(depthfile, name, haplogroup):
	with open(depthfile, 'r') as df:
		ddata = df.readlines()
	sc_depth = float(ddata[0].rstrip().split()[-1])
	if haplogroup == 'N/A':
		dfrow = [name, haplogroup, haplogroup, sc_depth]
	else:
		dfrow = [name, haplogroup[0], haplogroup, sc_depth]
	dfrow += [float(x) for x in ddata[2].rstrip().split()]
	return dfrow


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "Takes individual amplicon depth files and combines then into a pandas dataframe")

	parser.add_argument("infiles", nargs = "*", help = "individual files of normalized amplicon depths")
	parser.add_argument("-o", "--outfile", default = "./amplicon_depths_dataframe.pickle", help = "Output filename")
	parser.add_argument("-n", "--namelist", help = "File of male IDs of each file (default: filename of each file)")
	parser.add_argument("-g", "--haplogroups", help = "file of haplogroups of each individual")

	args = parser.parse_args()


	name_dict = {}
	if args.namelist:
		namesfile = open(args.namelist, 'r')
		names = namesfile.read().rstrip()
		names = names.split('\n')
		namesfile.close()
	
		if len(names) != len(args.infiles):
			print "ERROR: Number of input files and number of names do not match."
			quit()
		for idx, infile in enumerate(args.infiles):
			name_dict[infile] = names[idx]
	else:
		for infile in args.infiles:
			name_dict[infile] = infile.split('/')[-1]


	haplo_dict = {}
	if args.haplogroups:
		hapfile = open(args.haplogroups, 'r')
		haps = hapfile.read().rstrip()
		haps = {x.rstrip().split()[0]: x.rstrip().split()[1] for x in haps.split('\n')}
		hapfile.close()
		for infile in name_dict.keys():
			if name_dict[infile] in haps.keys():
				haplo_dict[infile] = haps[name_dict[infile]]
			else:
				haplo_dict[infile] = 'N/A'
	else:
		for infile in name_dict.keys():
			haplo_dict[infile] = 'N/A'
		
	column_names = ['Haplogroup', 'Sub-Haplogroup', 'Single-copy_depth', 'IR3', 'TSPY',
			'IR1', 'P8', 'P7', 'P6', 'P5', 'IR5', 'P4', 'DYZ19', 'RBMY', 'IR2', 'Blue',
			'Teal', 'Green', 'Red', 'DAZ', 'Gray', 'Yellow', 'Ctrl_reg_0', 'Ctrl_reg_1',
			'Ctrl_reg_2', 'Ctrl_reg_3', 'Ctrl_reg_4', 'Ctrl_reg_5']
			
	dfrows = []
	for infile in args.infiles:
		dfrows.append(depthfile_to_dfrow(infile, name_dict[infile], haplo_dict[infile]))
	dataframe = pd.DataFrame([x[1:] for x in dfrows], index = [x[0] for x in dfrows], columns=column_names)
	dataframe.to_pickle(args.outfile)
