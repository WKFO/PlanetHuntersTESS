import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator

cdpp_data = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_all_cdpp_data.txt', format='ascii',comment='#')


def cdpp(data, WF, n):

	'''
	
	function to output a new dataframe which includes the cdpp for each TIC (each lightcurve)

	input
	-----
	data (DF): data frame of aggregated data
	WF: which WF are we looking at (simply to keep track in the savign of files)
	n: either 05,10,20 - the range over which the cdpp was caluclated
	

	output 
	------

	DF with a cdpp column
	'''
	# where n in the time over which the cdpp was determined i.e. can be either 05,10,20

	obj_cdpp = cdpp_data['cdpp{}'.format(n)]
	obj_tic = cdpp_data['tic']

	cdpp = []

	obj_dict = {k: v for k, v in zip(obj_tic, obj_cdpp)}

	for i in data['TIC_ID']:

	    if int(i) in list(obj_tic):
	        cdpp.append(obj_dict[int(i)])

	    else:
	        cdpp.append(None)
	

	return pd.DataFrame(dict(CDPP=cdpp),index=data['TIC_ID'].index.values)


def df_extended_cdpp(data, WF,n):
	''''
	Function that merges the cdpp DF with the initial DF.

	'''
	cdpp_df = cdpp(data, WF, n)
	data_old = data.copy()
	data = data_old.merge(cdpp_df,how='outer', left_index=True, right_index=True)

	return data


def cdpp_stats_plot(data, WF, n, weight = None):

	'''
	Function that plots the cdpp for each TIC. 

	'''
	fig, ax = plt.subplots(figsize =(8,4))  # define the axes to plot

	data = df_extended_cdpp(data, WF,n)

	subject = ['planet', 'star', 'EB']
	color = ['green', 'orange', 'red']
	marker = ['s','o','v']


	for s,subj in enumerate(subject):
		# the convidence in the response when the majority was correct
		convidence = []
		cdpp = []
	

		if weight == None:
	
			for i,yes in enumerate(data['p_Yes']):
		
				if (data['subject'][i]) == subj:

					convidence.append(yes) # what percentage of people said 'yes' there is a planet
					cdpp.append(data['CDPP'][i]) # the cdpp of the star
			

		else:
	
			for i,yes in enumerate(data['p_Yes_weight']):
		
				if (data['subject'][i]) == subj:

					convidence.append(yes) # what percentage of people said 'yes' there is a planet
					cdpp.append(data['CDPP'][i]) # the cdpp of the star
			

		ax.plot(cdpp, convidence,marker=marker[s],markersize=4.5, ls='none',color=color[s],markeredgewidth=0.0, label = "{}".format(subject[s]))
	

	ax.set_xlabel('cdpp',fontsize=12)
	ax.set_ylabel('Convidence',fontsize=12)

	minorLocator = AutoMinorLocator()
	ax.xaxis.set_minor_locator(minorLocator)
	ax.tick_params(direction='in', which ='minor', colors='k',length = 5)
	minorLocator = AutoMinorLocator()
	ax.xaxis.set_minor_locator(minorLocator)
	ax.tick_params(direction='in', which ='minor', colors='k',length = 5)

	plt.legend(numpoints = 1)

	outpath = '/Users/Nora/Documents/research/TESS/ETE-6/output/results_figs/{}'.format(WF)
	plt.savefig(outpath + '/cdpp_stats_WF{}_{}{}.png'.format(WF,n,weight), format='png', dpi=80,bbox_inches='tight')






