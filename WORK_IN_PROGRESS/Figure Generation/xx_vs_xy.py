#!/usr/bin/python

'''
Make figure of XX vs. XY ampliconic coverage
'''


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
import glob
import numpy as np

# xx_vecfiles = glob.glob('/archive/page/2016.08.11-14522/solexa_page/lsteitz/XX_controls/???????/*vector')
# xy_vecfiles = glob.glob('/archive/page/2016.08.11-14522/solexa_page/lsteitz/XX_controls/XY_files/*vector')
xx_vecfiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/XX_controls_redo/???????_chr2_normed_amplicons')
xy_vecfiles = glob.glob('/lab/solexa_page/lsteitz/1000_Ys/XX_controls_redo/XY_???????_chr2_normed_amplicons')

ctrl_regs = 'IR3    TSPY    IR1     P8      P7      P6      P5      IR5     P4      DYZ19   RBMY    IR2     Blue    Teal    Green   Red     DAZ     Gray    Yellow  Ctrl_reg_0      Ctrl_reg_1      Ctrl_reg_2      Ctrl_reg_3      Ctrl_reg_4      Ctrl_reg_5'.split()
simple_ctrl_regs = 'IR3    IR1     P8      P7      P6      P5      IR5     P4     IR2     Blue    Teal    Green   Red      Gray    Yellow  Ctrl_reg_0      Ctrl_reg_1      Ctrl_reg_2      Ctrl_reg_5'.split()

xx_values = [[] for x in ctrl_regs]
xy_values = [[] for x in ctrl_regs]

for xxfile in xx_vecfiles:
    xx_data = open(xxfile, 'r')
    xx_vector = [float(x) for x in xx_data.readlines()[2].rstrip().split()]                 
    for amp_idx, amp_val in enumerate(xx_vector):
        xx_values[amp_idx].append(amp_val)
    xx_data.close()
    
for xyfile in xy_vecfiles:
    xy_data = open(xyfile, 'r')
    xy_vector = [float(x) for x in xy_data.readlines()[2].rstrip().split()]
    for amp_idx, amp_val in enumerate(xy_vector):
        xy_values[amp_idx].append(amp_val)
    xy_data.close()

plt.clf()
plt.figure(figsize=[20,5])
ax = plt.subplot(111)

for ampno, amp in enumerate(ctrl_regs):
    xys = plt.plot(np.arange(1.3, 1.8, 0.1) + ampno * 3,  xy_values[ampno], 'wo', markeredgecolor='black', label = 'XY Males' if ampno == 0 else '_')
    xxs = plt.plot(np.arange(0.8, 2.3, 0.1) + ampno * 3,  xx_values[ampno], 'kD', label = 'XX Females' if ampno == 0 else '_')

plt.xticks(np.arange(1.5, 75, 3), ctrl_regs[:19] + ['Ctrl %i' %(x) for x in range(1,7)])
plt.ylabel('Normalized coverage')

plt.tick_params(axis='x', direction='out', top='off')
plt.tick_params(axis='y', direction='out', right='off')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.axis([0,75,0,1.6])
plt.legend(numpoints = 1)
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/XX Controls/XX_vs_XY.png', dpi=300, bbox_inches='tight')
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/XX Controls/XX_vs_XY.svg', bbox_inches='tight')

plt.clf()
plt.figure(figsize=[7,1.75])
ax = plt.subplot(111)

xy_values = {ctrl_regs[x]:xy_values[x] for x in range(len(xy_values))}
xx_values = {ctrl_regs[x]:xx_values[x] for x in range(len(xx_values))}

for ampno, amp in enumerate(simple_ctrl_regs):
    xys = plt.plot(np.arange(1.3, 1.8, 0.1) + ampno * 3,  xy_values[amp], 'wo', markeredgecolor='black', label = 'XY Males' if ampno == 0 else '_')
    xxs = plt.plot(np.arange(0.8, 2.3, 0.1) + ampno * 3,  xx_values[amp], 'kD', label = 'XX Females' if ampno == 0 else '_')

plt.xticks(np.arange(1.5, 63, 3), simple_ctrl_regs[:15] + ['Ctrl %i' %(x) for x in range(1,7)], fontsize=7)
plt.ylabel('Normalized coverage', fontsize=10)

plt.tick_params(axis='x', direction='out', top='off')
plt.tick_params(axis='y', direction='out', right='off')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.axis([0,63,0,1])
lgd = plt.legend(numpoints = 1, ncol=2, fontsize=10)
lgd.get_frame().set_linewidth(0.5)
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/XX Controls/XX_vs_XY_simple_regs.png', dpi=300, bbox_inches='tight')
plt.savefig('/lab/Page_lab-users/lsteitz/1000_Genomes_Male_Figures/XX Controls/XX_vs_XY_simple_regs.svg', bbox_inches='tight')
