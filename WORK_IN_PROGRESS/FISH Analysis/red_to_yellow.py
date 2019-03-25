#!/usr/bin/python

'''
Takes FISH .tif files and turns the red dots yellow
'''

import sys
import numpy as np
from PIL import Image

#filter_value = 70

for imagefile in sys.argv[1:]:
	if not imagefile.endswith('.tif'):
		print 'ERROR: Image does not have .tif suffix.'
		quit()
	image = Image.open(imagefile)
	r, g, b = image.split()
	rar = np.array(r)
	gar = np.array(g)
	bar = np.array(b)

	
	image_len = len(rar)
	image_wid = len(rar[0])
	
	for pix_l in xrange(image_len):
		for pix_w in xrange(image_wid):
			if gar[pix_l][pix_w] < rar[pix_l][pix_w]:
				gar[pix_l][pix_w] = rar[pix_l][pix_w]


	rf = Image.fromarray(rar)
	gf = Image.fromarray(gar)
	bf = Image.fromarray(bar)

	filtered_image = Image.merge('RGB', (rf, gf, bf))
	
	imagename = imagefile.split('/')[-1]
	outname = '/Users/lsteitz/Desktop/' + imagename[:-4] + '_yellowed.tif'
	filtered_image.save(outname)
