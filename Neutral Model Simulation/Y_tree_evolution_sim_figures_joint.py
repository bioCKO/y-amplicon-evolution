#!/usr/bin/python

'''
Makes single figure of one-parameter model across all parameters
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import glob
import numpy as np

real_muts = 204
real_fitch = 116
#real_muts = 183
#real_fitch = 103

color_norm = matplotlib.colors.LogNorm(vmin=.000001, vmax=.005)
img = plt.imshow(np.array([[0,1]]), norm=color_norm, cmap="inferno")
img.set_visible(False)

plt.clf()
plt.figure(figsize=[10,6])
ax = plt.subplot(111)
model_mean_muts = []
model_mean_fitch = []

plt.plot(real_muts, real_fitch, 'D', mfc='#6B8392', mec='#596C78', mew=0.7, ms=10, zorder=10)


infiles = ['/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du5_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du3_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du2_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du1_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du05_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du03_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du02_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du01_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du005_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du003_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du002_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du001_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du0005_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du0003_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du0002_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du0001_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du5e-05_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du3e-05_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du2e-05_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du1e-05_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du5e-06_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du3e-06_de0_dur0_10000_sims.txt', 
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du2e-06_de0_dur0_10000_sims.txt',
'/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/final_sim_stats_no_gen_estimate/du1e-06_de0_dur0_10000_sims.txt'] 


for input_file in infiles[8:]:
	
	duval = input_file.split('/')[-1].split('_')[0][2:]
	if not (duval == '0' or 'e' in duval):
		duval = '.' + duval
	duval = float(duval)	
	if duval == 0.:
		continue
	
	sim_muts = []
	sim_fitch = []
	datafile = open(input_file, 'r')
	for line in datafile:
		data = map(int, line.rstrip().split())
		sim_muts.append(data[0])
		sim_fitch.append(data[1])
	datafile.close()

	plt.plot(sim_muts, sim_fitch, '.', color=matplotlib.cm.inferno(color_norm(duval)), zorder=1, alpha = 0.75, mew=0.1)

	model_mean_muts.append(np.mean(sim_muts))
	model_mean_fitch.append(np.mean(sim_fitch))

#plt.plot(model_mean_muts, model_mean_fitch, 'o', color='black', zorder=2)
plt.plot([], [], 'D', color='black', ms=10, label='Real data')
plt.plot([], [], '.', color='black', label='Simulated data')
plt.xlabel('Number of mutants', size=14)
plt.ylabel('Number of mutation events', size=14)

cb = plt.colorbar(img, label='Mutation rate', orientation='vertical')
cb.set_ticks([.000001, .00001, .0001, .001])


#matplotlib.colorbar.ColorbarBase(cmap="inferno", norm=color_norm, orientation="vertical")

#plt.title('Single-parameter model')
lgd = plt.legend(numpoints=1, loc='center', bbox_to_anchor=(0.5,-0.17), ncol=2)
lgd.get_frame().set_linewidth(0.)
#plt.colorbar()
plt.axis([0,1200,0,130])
plt.tick_params(axis='x', direction='out', top='off')
plt.tick_params(axis='y', direction='out', right='off')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
#plt.tight_layout()

outname = '/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Phylogenetic_Trees/Simulations/one_par_model_joint_10000_sims_no_gen_estimate'
plt.savefig('%s.png' %(outname), bbox_extra_artists=[lgd], bbox_inches='tight', dpi=300)
plt.savefig('%s.svg' %(outname), bbox_extra_artists=[lgd], bbox_inches='tight')