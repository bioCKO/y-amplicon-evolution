#!/usr/bin/python

from make_chrY_vector import bed_to_depthvec
#get regs

rfile = open('/lab/Page_lab-users/lsteitz/Y_repeats.txt', 'r')

y_reps = {}
for line in rfile:
	if line.startswith('#'):
		continue
	data = line.rstrip().split('\t')
	if not data[5] in y_reps:
		y_reps[data[5]] = []
	y_reps[data[5]].append(data)

rfile.close()

#source_bed = 'chrY_to_all_100bp_cut14.bedgraph'
#source_bed = 'chrY_to_all_100bp_cut0.bedgraph'
#source_bed = 'chrY_to_all_nors_100bp_cut0.bedgraph'
#source_bed = 'chrY_to_all_norsplusten_35bp_cut14.bedgraph'

source_beds = ['chrY_to_all_100bp_cut14.bedgraph', 
			   'chrY_to_all_100bp_cut11.bedgraph', 
			   'chrY_to_all_100bp_cut0.bedgraph', 
			   'chrY_to_all_50bp_cut14.bedgraph', 
			   'chrY_to_all_50bp_cut11.bedgraph', 
			   'chrY_to_all_50bp_cut0.bedgraph', 
			   'chrY_to_all_35bp_cut14.bedgraph', 
			   'chrY_to_all_35bp_cut11.bedgraph', 
			   'chrY_to_all_35bp_cut0.bedgraph']

print '\t' + '\t'.join(['_'.join(x.split('.')[0].split('_')[3:5]) for x in source_beds])

y_rep_pcts = {x:[] for x in y_reps}

for source_bed in source_beds:
	winlen = int(source_bed.split('_')[3][:-2])
	bases = bed_to_depthvec(source_bed, 'chrY')

	for reg in y_reps:
		#print reg,
		total_good_bases = 0
		total_bases = 0
		for subreg in y_reps[reg]:
			good_bases = 0
			#for base in bases[int(subreg[1]):int(subreg[2])]:
				#if 0 < base <= int(subreg[8]):
			for base in xrange(int(subreg[1]),int(subreg[2])):
				if all([0 < x <= int(subreg[8]) for x in bases[base - winlen + 1:base + 1]]):
					good_bases += 1
			#print '\t\t', good_bases, int(subreg[7]), good_bases / float(subreg[7])
			total_good_bases += good_bases
			total_bases += int(subreg[7])
		y_rep_pcts[reg].append(total_good_bases / float(total_bases))
		#print '\t', total_good_bases, total_bases, total_good_bases / float(total_bases)

for reg in y_reps:
	print '%s\t%s' %(reg, '\t'.join(map(str, y_rep_pcts[reg])))