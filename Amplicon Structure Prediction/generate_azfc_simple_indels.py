#!/usr/bin/python

azfc_ref = 'bttb1grreby5grrg5ybe'

amps = '15btgrey'

cns = set([])

for dstart in xrange(0, len(azfc_ref)):
	for dend in xrange(dstart + 1, len(azfc_ref)+1):
		azfc_del = azfc_ref[:dstart] + azfc_ref[dend:]
		azfc_dup = azfc_ref[:dend] + azfc_ref[dstart:]
		del_counts = [azfc_del.count(x) for x in amps]
		dup_counts = [azfc_dup.count(x) for x in amps]
		del_counts = tuple(del_counts)
		dup_counts = tuple(dup_counts)
		cns.add(del_counts)
		cns.add(dup_counts)

		
		
		
