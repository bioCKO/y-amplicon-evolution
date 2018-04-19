#!/usr/bin/python

'''
Get coverage stats of regions in a bedfile, given a coverage file
'''

import argparse
import numpy
import scipy.stats

parser = argparse.ArgumentParser(description = "Get coverage stats of regions in a bedfile, given a coverage file")
parser.add_argument('coverage', help = "Coverage file")
parser.add_argument('regions', help = "BED file of regions")
parser.add_argument('-e', '--edge', help = "Ignore first and last X bases of regions, to account for edge effects")

args = parser.parse_args()

with open(