#!/usr/bin/python

import sys
import bed_to_vector as btv

ctrl = btv.bed_to_depthvec(sys.argv[1], 'chrY')
cnvec = btv.bed_to_depthvec(sys.argv[2], 'chrY')

startbase = 3000000
endbase = 9000000
bases = range(3000000,9000000)+range(12000000,26000000)

locs = 0
ctrlcn = 0
cn = 0
difbases = 0
onebases =  0
onediffs = 0

for base in bases:
	if ctrl[base] > 0:
		locs += 1
		ctrlcn += ctrl[base]
		cn += cnvec[base]
		if ctrl[base] != cnvec[base]:
			difbases += 1
	if ctrl[base] == 1:
		onebases += 1
		if cnvec[base] != 1:
			onediffs += 1

print 'Mappable locations: %i' %(locs)
print 'Control mean cn: %.4f' %(float(ctrlcn)/locs)
print 'Mean cn: %.4f' %(float(cn)/locs)
print 'Percent bases with difference: %.4f' %(float(difbases)/locs)
print 'Percent 1-cn bases in ctrl with >1 cn: %.4f' %(float(onediffs)/onebases)
