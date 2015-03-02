from __future__ import print_function
from __future__ import division

import sys
import os
import numpy as np

#import logging
#import pdb

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

#import pandas as pd
import vcf

# Insert the src/python/directory at front of the path
rvddir = os.path.join('../../src/python/rvd27')
sys.path.insert(0, rvddir)
#import rvd27

def main():
	vcfFile = 'Varscan2_somatic_HCC1187_expanded.vcf'
	vcf_reader = vcf.Reader(open(vcfFile, 'r'))
	varscan_position =[record.POS for record in vcf_reader]

        muTect_position = [154736692, 154743899,154749704,154753635, 154754371,154758813, 154760439,154766700,
            154766732,154777014,154777118,154780960,154784369,154787357,154787367,154794220]
        
	RVD2_somatic_position = [154749704, 154753635,154754371,154758813,
                                 154760439,154766732,154766832,154777118,154787612,154788798]
	
	RVD2_germline_position = [154743899,154749704, 154753635,154754371,154758813,
                                   154780960,154784369, 154787357,154787367, 154787612,154788798]

	Strelka_position = [154760439]

	position = list(set(varscan_position+muTect_position+
		RVD2_somatic_position+RVD2_germline_position+Strelka_position))

	position = sorted(position)

	# plot varscan
	varscan_idx = [position.index(pos)+1 for pos in varscan_position]

	muTect_idx = [position.index(pos)+1 for pos in muTect_position]

	RVD2_somatic_idx = [position.index(pos)+1 for pos in RVD2_somatic_position]

	RVD2_germline_idx = [position.index(pos)+1 for pos in RVD2_germline_position]

	Strelka_idx = [position.index(pos)+1 for pos in Strelka_position]

	fig = plt.figure(figsize=(12,4))
	gs = gridspec.GridSpec(1, 2, width_ratios=[5.2, 1.5]) 
	
	# ax0
	ax0 = plt.subplot(gs[0])
	ax0. plot(varscan_idx, np.ones_like(varscan_idx),'r*')
	ax0. plot(muTect_idx, 2*np.ones_like(muTect_idx),'md')
	ax0. plot(RVD2_somatic_idx, 3*np.ones_like(RVD2_somatic_idx),'yo')
	ax0. plot(RVD2_germline_idx, 4*np.ones_like(RVD2_germline_idx),'bx')
	ax0. plot(Strelka_idx, 5*np.ones_like(Strelka_idx),'gv')

	Idx = [1,8,12,15,17,24,26,39,41,55,60,71,73,75,77,84]

	ax0.text(Idx[0]-2.9, 5.3,'Idx:%d'%Idx[0],fontsize=9, backgroundcolor = 'w')
	for i in xrange(len(Idx)-1):
		ax0.text(Idx[i+1]-0.7, 5.3,'%d'%Idx[i+1],fontsize=9, backgroundcolor = 'w')

	ax0.set_ylim([0,6])
	ax0.set_xticks(np.arange(5, len(position)+5, 5.0))
	label = ['','VarScan2\nSomatic ','muTect ', 'RVD2\nSomatic ',
                 'RVD2\nGermline ', 'Strelka ','']
	ax0.set_yticklabels(label)
	ax0.set_xlabel('Position Index')

	ml = MultipleLocator(1)
	ax0.xaxis.set_minor_locator(ml)

	ax0.grid(True,which='both')

        # ax1
	ax1 = plt.subplot(gs[1],sharey=ax0)
	xval = [len(varscan_idx), len(muTect_idx),
	 len(RVD2_somatic_idx), len(RVD2_germline_idx), len(Strelka_idx)]
	ypos = [1,2,3,4,5]

	ax1.barh(ypos[0],xval[0], color='r', align='center',height=0.1) 
	ax1.text(xval[0]+2,ypos[0]-0.05,'%d' %xval[0])
	ax1.barh(ypos[1],xval[1], color='m', align='center',height=0.1)
	ax1.text(xval[1]+2,ypos[1]-0.05,'%d' %xval[1])
	ax1.barh(ypos[2],xval[2], color='y', align='center',height=0.1)
	ax1.text(xval[2]+2,ypos[2]-0.05,'%d' %xval[2])
	ax1.barh(ypos[3],xval[3], color='b', align='center',height=0.1)
	ax1.text(xval[3]+2,ypos[3]-0.05,'%d' %xval[3])
	ax1.barh(ypos[4],xval[4], color='g', align='center',height=0.1)
	ax1.text(xval[4]+2,ypos[4]-0.05,'%d' %xval[4])
	ax1.set_xlim([0,95])
	# make these tick labels invisible
	plt.setp( ax1.get_yticklabels(), visible=False)
	ax1.set_xlabel('Total')

	plt.tight_layout()
	# plt.show()
	plt.savefig('HCC1187_expanded_dotplot.pdf')

	## write positions to a file
	File = 'expanded_position_lookup_chart.txt'
	f = open(File,'w')
	print('Index\tActual Position', file=f)
	for i in xrange(len(position)):
		print('%d\tchr7:%d'%(i+1, position[i]), file=f)

	f.close()


if __name__ == "__main__":
	main()
