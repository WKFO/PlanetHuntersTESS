import numpy as np
import pandas as pc

import matplotlib.pyplot as plt


def SNR_stats(data, WF, weight = None):

	# the convidence in the response when the majority was correct and SNR of planet
	convidence_correct = []
	SNR_correct = []

	# the convidence in the response when the majority was wrong and SNR of planet
	convidence_false = []
	SNR_false = []

	# only use the identifiations of the planets (don't have SNR data on the rest)
	#only_planets = data.subject_type == "planet"
	#data = data[only_planets] # mask the rest
	
	if weight == None:

		for i,yes in enumerate(data['p_Yes']):
	
			if (data['subject'][i]) == 'planet':
				# there is a planet
				if yes > data['p_No'][i]:
					convidence_correct.append(yes)
					SNR_correct.append(data['SNR'][i])
		
				if yes < data['p_No'][i]:
					convidence_false.append(yes)
					SNR_false.append(data['SNR'][i])
	else:

		for i,yes in enumerate(data['p_Yes_weight']):
	
			if (data['subject'][i]) == 'planet':
				# there is a planet
				if yes > data['p_No_weight'][i]:
					convidence_correct.append(yes)
					SNR_correct.append(data['SNR'][i])
		
				if yes < data['p_No_weight'][i]:
					convidence_false.append(yes)
					SNR_false.append(data['SNR'][i])


	fig, ax = plt.subplots(figsize =(10,4))  # define the axes to plot

	ax.plot(SNR_correct, convidence_correct,marker='o',markersize=4.5, ls='none',color='green',markeredgewidth=0.0)
	ax.plot(SNR_false, convidence_false,marker='x',markersize=4.5, ls='none',color='red',markeredgewidth=1.0)

	ax.set_xlabel('SNR',fontsize=12)
	ax.set_ylabel('Convidence',fontsize=12)

	minorLocator = AutoMinorLocator()
	ax.xaxis.set_minor_locator(minorLocator)
	ax.tick_params(direction='in', which ='minor', colors='k',length = 5)
	minorLocator = AutoMinorLocator()
	ax.xaxis.set_minor_locator(minorLocator)
	ax.tick_params(direction='in', which ='minor', colors='k',length = 5)

	plt.savefig('results_figs/SNR_stats_WF{}{}.png'.format(WF,weight), format='png', dpi=80,bbox_inches='tight')

#SNR_stats(1, weight = "w")
#SNR_stats(2, weight = "w")


print("Getting ground truth information for SNR values...")

def SNR(data, WF):


	pl_snr = pl['col16']
	pl_tic = pl['col1']

	snr = []

	pl_dict = {k: v for k, v in zip(pl_tic, pl_snr)}
	
	for i in data['TIC_ID']:

	    if int(i) in list(pl_tic):
	        snr.append(pl_dict[int(i)])
	    else:
	        snr.append(None)
	
	return pd.DataFrame(dict(SNR=snr),index=data['TIC_ID'].index.values)
