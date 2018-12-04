import numpy as np
import pandas as pc

import matplotlib.pyplot as plt

outpath = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/figures/'

def classification_count(data, weight = None):
	'''
	Count how many people said there were planets when there were/weren't etc. 

	input
	-----
	data (DF): data frame of aggregated data
	WF: which WF are we looking at (simply to keep track in the savign of files)
	weight: weighting on or off (off by default)

	output
	-----
	true_planet, true_noplanets, false_planets, false_noplanets, true_eb, false_eb

	'''
	true_transit = 0
	false_transit = 0 # there are planets but they said no

	for i,subject in enumerate(data['subject']):

		# there is a planet
		if weight == None:

			if subject == 'planet' and (data['p_Yes'][i] > data['p_No'][i]):
				true_transit +=1
			elif subject == 'planet' and (data['p_Yes'][i] < data['p_No'][i]):
				false_transit +=1	
	
	
			#there is an EB
			elif subject == 'EB' and (data['p_Yes'][i] > data['p_No'][i]):
				true_transit +=1	
			elif subject == 'EB' and (data['p_Yes'][i] < data['p_No'][i]):
				false_transit +=1
	
		# there is a planet
		else:

			if subject == 'planet' and (data['p_Yes_weight'][i] > data['p_No_weight'][i]):
				true_transit +=1
			elif subject == 'planet' and (data['p_Yes_weight'][i] < data['p_No_weight'][i]):
				false_transit +=1	
	
			#there is an EB
			elif subject == 'EB' and (data['p_Yes_weight'][i] > data['p_No_weight'][i]):
				true_transit +=1	
			elif subject == 'EB' and (data['p_Yes_weight'][i] < data['p_No_weight'][i]):
				false_transit+=1	
	
	return true_transit, false_transit  # we don't have informatin on them marking things when there's nothing there as only give them sims with planets.



def pie_identification_plot(data):

	'''
	
	Pie chart of how many people identified planets correctly, incorrectly. 

	input
	-----
	data (DF): data frame of aggregated data
	WF: which WF are we looking at (simply to keep track in the savign of files)
	weight: weighting on or off (off by default)
	
	output
	------
	3 pie chart. 

	'''

	#num_pl = data.subject.value_counts()['planet']
	#num_eb = data.subject.value_counts()['EB']
	#num_star = data.subject.value_counts()['star']
	
	true_transit, false_transit = classification_count(data, weight= None)
	true_transit_w, false_transit_w = classification_count(data, weight = 'w')

	# Pie chart, where the slices will be ordered and plotted counter-clockwise:
	#labels = 'Correct', 'Wrong'
	colors ='limegreen', 'orangered'
	sizes_p = [true_transit, false_transit]


	f, (ax1, ax2) = plt.subplots(1,2, figsize =(10,5))  # define the axes to plot
	
	ax1.set_title("Transit")
	ax1.pie(sizes_p, colors =colors, autopct='%1.1f%%',
	        shadow=True, startangle=90, wedgeprops = {'linewidth': 0})
	ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

	sizes_p_w = [true_transit_w, false_transit_w]

	ax2.set_title("Transit w")
	ax2.pie(sizes_p_w, colors =colors, autopct='%1.1f%%',
	        shadow=True, startangle=90, wedgeprops = {'linewidth': 0})
	ax2.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	

	
	plt.savefig(outpath + '/pie_chart_sims.png', format='png', dpi=80)









