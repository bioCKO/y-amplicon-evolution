#!/usr/bin/python

'''
Takes FISH .tif files and filters out low-intensity pixels in the RGB bands of the images
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
	
	# r_filter = np.mean(rar.flatten()), np.std(rar.flatten())
# 	g_filter = np.mean(gar.flatten()), np.std(gar.flatten())
# 	b_filter = np.mean(bar.flatten()), np.std(bar.flatten())
# 	
# 	print r_filter, g_filter, b_filter
	
	r_filter = 2.2 * np.mean(rar.flatten())
	g_filter = 2 * np.mean(gar.flatten())
	b_filter = 2 * np.mean(bar.flatten())

	
	image_len = len(rar)
	image_wid = len(rar[0])
	
	for pix_l in xrange(image_len):
		for pix_w in xrange(image_wid):
			if rar[pix_l][pix_w] < r_filter:
				rar[pix_l][pix_w] = 0
			if gar[pix_l][pix_w] < g_filter:
				gar[pix_l][pix_w] = 0
			if bar[pix_l][pix_w] < b_filter:
				bar[pix_l][pix_w] = 0

	rf = Image.fromarray(rar)
	gf = Image.fromarray(gar)
	bf = Image.fromarray(bar)

	filtered_image = Image.merge('RGB', (rf, gf, bf))
	
	outname = imagefile[:-4] + '_filtered.tif'
	filtered_image.save(outname)
