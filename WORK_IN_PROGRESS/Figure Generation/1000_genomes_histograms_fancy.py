#!/usr/bin/python

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker

vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_amplicons.pickle')
#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls_badctrlfiltered.pickle')
#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_gccorrected_with_controls.pickle')
#vectors = pd.read_pickle('/lab/solexa_page/lsteitz/1000_Ys/1000_Genomes_simple_Y_repeats_vectors_10bpbuffer.pickle')

# import matplotlib.font_manager
# print matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')


#sns.plotting_context('poster')
#sns.set_color_codes()
max_depth = 2
'''
for color in [['Blue', 'b', 4.], ['Red', 'r', 4.], ['Yellow', 'y', 2.], ['Green', 'g', 3.], ['P8', '#70D6D0', 2.]]: #[['IR1', '#30C4BB', 2.]]:#
	hist_data = plt.hist(vectors[color[0]], bins=np.arange(0,max_depth,.01))
	max_val = max(hist_data[0])
	plt.clf()
	plt.figure(figsize=[4,3])
	#Big text sizes for full size figure: x-axis test: 18; y-axis text: 18; ticks: 13; title: 24
	ax = sns.distplot(vectors[color[0]], bins=np.arange(0,max_depth,.01), hist=True, kde=False, color=color[1])
	plt.xticks(np.arange(1/(color[2]*2), max_depth, 1/color[2]), [])
	plt.yticks([])
	plt.text(max_depth/2., max_val/-9, 'Normalized depth', size=12,
			 horizontalalignment='center', verticalalignment='center')

	plt.text(-.205, max_val*.55, 'Number of men', size=12,
			 horizontalalignment='center', verticalalignment='center', rotation='vertical')
	for xtick in np.arange(0, max_depth + 0.1, .25):
		plt.text(xtick, max_val/-78, str(xtick), size=9,
		horizontalalignment='center', verticalalignment='top')
	for ytick in np.arange(0, max_val * 1.1, 25):
		plt.text(-0.02, ytick, str(int(ytick)), size=9,
		horizontalalignment='right', verticalalignment='center')
	ax.axis([0,max_depth,0,max_val * 1.1])
	plt.text(0.019, .85 * max_val, '0 copies', 
			 horizontalalignment='left', 
			 verticalalignment='bottom', 
			 rotation='vertical',
			 color='#606060', 
			 size=9)
	plt.text(1/(color[2]*2) + 0.019, .85 * max_val, '1 copy', 
			 horizontalalignment='left', 
			 verticalalignment='bottom', 
			 rotation='vertical',
			 color='#606060', 
			 size=9)		 
	for cn_index, copy_number in enumerate(np.arange(3/(color[2]*2), max_depth, 1/color[2])):
		plt.text(copy_number + 0.019, .85 * max_val, '%i copies' %(cn_index + 2), 
				 horizontalalignment='left', 
				 verticalalignment='bottom', 
				 rotation='vertical',
				 color='#606060', 
				 size=9)
	plt.xlabel('')
	plt.title('%s Amplicon' %(color[0]), size=14)
	#plt.plot([.98263,.98263],[0,1000],'black')
	
	
	#import matplotlib.mlab as mlab
	
	
	#plt.plot(np.linspace(0,2,500), mlab.normpdf(np.linspace(0,2,500), np.mean(vectors[color[0]]), np.std(vectors[color[0]])))
	
	
	
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Simple Repeat Histograms Fancy/%s_Copy_Number_Histogram_bigtext.svg' %(color[0]))#, dpi=300)
	

'''
for ctrl_region in xrange(6):
	hist_data = plt.hist(vectors['Ctrl_reg_%i' %(ctrl_region)], bins=np.arange(0,max_depth,.01))
	max_val = max(hist_data[0])
	plt.clf()
	plt.figure(figsize=[2, 1.5])
	ax = sns.distplot(vectors['Ctrl_reg_%i' %(ctrl_region)], bins=np.arange(0,max_depth,.01), hist=True, kde=False)#, color=color[1])
	plt.xticks([])
	plt.yticks([])
	plt.text(max_depth/2., max_val/-6, 'Normalized depth', size=8,
			 horizontalalignment='center', verticalalignment='center')

	plt.text(-.3, max_val*.55, 'Number of men', size=8,
			 horizontalalignment='center', verticalalignment='center', rotation='vertical')
	for xtick in np.arange(0, max_depth + 0.1, .5):
		plt.text(xtick, max_val/-78, str(xtick), size=6,
		horizontalalignment='center', verticalalignment='top')
	for ytick in np.arange(0, max_val * 1.1, 25):
		plt.text(-0.02, ytick, str(int(ytick)), size=6,
		horizontalalignment='right', verticalalignment='center')
	ax.axis([0,max_depth,0,max_val * 1.1])
	plt.xlabel('')
	plt.title('Control Region %i' %(ctrl_region + 1), size=10)
	plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/Simple Repeat Histograms Fancy/Control_%i_Copy_Number_Histogram_forpaper.svg' %(ctrl_region))#, dpi=300)
