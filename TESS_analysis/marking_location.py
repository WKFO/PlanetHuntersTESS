
import numpy as np
import pandas as pd 
from glob import glob
import astropy.io.fits as pf

from astropy.table import Table
import matplotlib.pyplot as plt
import ast


outpath = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/figures/'



def string_to_list(vals):
	'''
	The input values are in a string which look like a list. Convert this into something that python recognises as a list. 

	Input:
	----
	vals (str): string that looks like a list

	oputput:
	------
	vals (list)
	'''

	# split the values up into indivdual parts, each of which is still a string
	try:
		vals = ast.literal_eval(vals)
		vals = [float(v) for v in vals]
	except:
		vals = [np.nan, np.nan]
	
	return vals


def pixel_to_time(xvals, tmin, tmax, xmin, xmax):
	"""
	Convert the x axis from pixels into hours. 

	The number in this function will change depending on what the conversio between in the images and the axes are.
	

	input:
	-----
	xvals (list): list of the x values ( in this case for a given TIC) in pixels
	
	output:
	------
	x_vals (list) x values in days.
	"""
	# calculate the conversion factor

	x_time = []

	for x in xvals:

		x_t = ((x - xmin) / (1193)) * (tmax - tmin)
		
		x_time.append(x_t)
		
	return x_time


def time_to_pixel(time, tmin, tmax):
	"""
	Convert the x axis from pixels into hours. 

	The number in this function will change depending on what the conversio between in the images and the axes are.
	

	input:
	-----
	xvals (list): list of the x values ( in this case for a given TIC) in pixels
	
	output:
	------
	x_vals (list) x values in days.
	"""
	x_pix = []

	for t in time:

		x_pixel = ((t - tmin) / (tmax - tmin)) * 1193 # not sure why this is 

		x_pix.append(x_pixel)


	return x_pix



pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data.txt', format='ascii',comment='#')
eb = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_eb_data.txt', format='ascii',comment='#')



#transit_peaks_hours(data_WF)
def between(l1,low,high):
    l2 = []
    for i in l1:
        if(i > low and i < high):
            l2.append(i)
    return l2


# ----------------------------------------

# plot the original image and the distribution of markings

def marking_plot(data, sort = 'p_Yes', dist = 'gauss', number = 20):

	'''
	This function overplots the locations of the marking ontop of the LC png. This will allow us to locate where people
	have marked transits and quickly verify wherthere these are true transits. 

	Will want to plot these for LC with a certain 'likelihood' of a transit. 


	input 
	-----
	data : the data frame with all the data (output from the extraction code)
	sort : what to sort the data by in descending order. choices: p_Yes, p_Yes_weight,x_peak_max,x_peak_max_box
	dist: (short for distribution), can be either 'gauss'for the gaussian distribution or 'box' for the 'box' distribution
	i.e. how each transit is modelled. Both are weighted. 
	number: how many plots to make. By default it will plot the top 20 transit candidates.

	'''
	# sort the data frame to have the most likely transits at the top
	# can sort eithe rby yes/no, peak box or by peak gauss

	# REMEMBER to keep updating this list
	mast_planets = Table.read('/Users/Nora/Documents/research/TESS/planethunters/Rel01/MAST_planets/MAST_planets.txt', format='ascii',comment='#')


	data = data.sort_values(by = sort, ascending = False)

	# select the top 'number' of candidates

	data = data.head(n=number)

	for rown_n, row in data.iterrows(): # do for each row indivisually



		# CHANGE TO FALSEEEE !!!!!!!!!!!!! OTHERWISE ONLY LOOKING AT SIMS DATA THIS IS VERY IMPORTANT AAAH



		if row['subject_type'] == True: # only plot if not already a known planet. Change to assess nby 'sim'

			image_ID = row['candidate']
			
			img = plt.imread("/Users/Nora/Documents/research/TESS/planethunters/Rel01/lightcurves/data/{}".format(image_ID))

			fig, ax = plt.subplots()
			ax.imshow(img, extent=[0, 1193, 0, 500])
	
			xmin = row['min_time']
			xmax = row['max_time']
			
			xminpix = row['min_pix']
			xmaxpix = row['max_pix']


			if dist == 'gauss':
				pdf_peaks = row['pdf']       
				transit_pix = row['x_eval']  


			elif dist == 'box':
				pdf_peaks = row['pdf_box']       
				transit_pix = row['x_eval_box']  				
			

			transit_pix = string_to_list(transit_pix)
			pdf_peaks = string_to_list(pdf_peaks)
			
			pdf_peaks = [i * (400/np.max(pdf_peaks)) for i in pdf_peaks]

			if row['TIC_ID'] in mast_planets['tic_id']:
				plt.title("TIC ID: {} has a MAST confirmed planet. {} Rank: {}".format(row['TIC_ID'],dist,rown_n + 1))
			else:
				plt.title("TIC ID: {} has no known planets yet. {} Rank: {}".format(row['TIC_ID'],dist, rown_n + 1))

			plt.plot(transit_pix, pdf_peaks, '--', linewidth=1, color='firebrick', alpha = 0.9)
			
			plt.savefig(outpath + '/marked_transit_distribution_{}_{}_{}.png'.format(dist,rown_n + 1,row['TIC_ID']), format='png', dpi=300,bbox_inches='tight')



def marking_method_comp(data, sort = 'p_Yes', dist = 'all', number = 20):

	'''
	This function plots the returned most likely transits. 
	Compared box fitting, gaussian fitting and db scan to find locations. 
	Plot of them on the original image. 
	'''

	# Import published list of MAST found planets

	mast_planets = Table.read('/Users/Nora/Documents/research/TESS/planethunters/Rel01/MAST_planets/MAST_planets.txt', format='ascii',comment='#')

	# sort the data - we want to try this for the most likely candidates.
	data = data.sort_values(by = sort, ascending = False)

	# select the top n candidates which are the most likely to have a transit.
	data = data.head(n=number)


	for rown_n, row in data.iterrows(): # do for each row indivisually



		# CHANGE TO FALSEEEE !!!!!!!!!!!!! OTHERWISE ONLY LOOKING AT SIMS DATA THIS IS VERY IMPORTANT AAAH




		if row['subject_type'] == True: # only plot if not a sim. Change to assess by 'sim' for real data.

			image_ID = row['candidate']
			
			img = plt.imread("/Users/Nora/Documents/research/TESS/planethunters/Rel01/lightcurves/data/{}".format(image_ID))

			fig, ax = plt.subplots()
			ax.imshow(img, extent=[0, 1193, 0, 500])

			xmin = row['min_time']
			xmax = row['max_time']
			
			xminpix = row['min_pix']
			xmaxpix = row['max_pix']
		
			#pdf_peaks = row['pdf']       
			#transit_pix = row['x_eval']  

			pdf_peaks = string_to_list(row['pdf_peaks'])
			x_peaks = string_to_list(row['x_peaks'])
			pdf_peaks_box = string_to_list(row['pdf_peaks_box'])
			x_peaks_box = string_to_list(row['x_peaks_box'])
			db_peak = string_to_list(row['db_peak'])
			db_count = string_to_list(row['db_count'])
			db_count_weighted = string_to_list(row['db_count_weighted'])

			# scale them to make the appear better on the plot, height represents strength of the vote
			pdf_peaks = [i * (400/np.max(pdf_peaks)) for i in pdf_peaks]
			pdf_peaks_box = [i * (400/np.max(pdf_peaks_box)) for i in pdf_peaks_box]
			db_count = [i * (400/np.max(db_count)) for i in db_count]
			db_count_weighted = [i * (400/np.max(db_count_weighted)) for i in db_count_weighted]
			

			if row['TIC_ID'] in mast_planets['tic_id']:
				plt.title("TIC ID: {} has a MAST confirmed planet. {} Rank:".format(row['TIC_ID'],dist,rown_n + 1))
			else:
				plt.title("TIC ID: {} has no known planets yet. {} Rank: {}".format(row['TIC_ID'],dist, rown_n + 1))


			if dist == 'all':
				plt.plot(x_peaks, pdf_peaks, marker = 'o', markersize = 12, ls = 'None', color='firebrick', label = 'gaussian')
				plt.plot(x_peaks_box, pdf_peaks_box, marker = 'x', markersize = 12,ls = 'None',  color='c', label = 'box')
				plt.plot(db_peak, db_count, marker = 's', markersize = 12, ls = 'None', color='m', label = 'DB')
				plt.plot(db_peak, db_count_weighted, marker = '<', markersize = 12, ls = 'None', color='yellow', label = 'DB weighted')

			elif dist == 'gauss':
				plt.plot(x_peaks, pdf_peaks, marker = 'o', markersize = 12, ls = 'None', color='firebrick', label = 'gaussian')
			elif dist == 'box':
				plt.plot(x_peaks_box, pdf_peaks_box, marker = 'x', markersize = 12,ls = 'None',  color='c', label = 'box')
			elif dist == 'DB':
				plt.plot(db_peak, db_count, marker = 's', markersize = 12, ls = 'None', color='m', label = 'DB')
				plt.plot(db_peak, db_count_weighted, marker = '<', markersize = 12, ls = 'None', color='yellow', label = 'DB weighted')


			plt.legend(numpoints=1, fontsize = 10, bbox_to_anchor=(1, -0.1))
			
			plt.savefig(outpath + '/marked_transit_peaks_{}_{}_{}.png'.format(dist,rown_n + 1,row['TIC_ID']), format='png', dpi=300,bbox_inches='tight')



# --------------
# May use later 


def peak_location_column(data, WF):

	'''
	# compare the marked transit to the injected transit
	# model the injected transit as a box, width = duration of transit (column 9 of table)
	# model each marked transit as a box too and see whether they overlap and later see if the centre of the marked box overlaps?

	input
	-----
	data (DF): data frma of all the data
	tolerance (float): the +/- tolerance in days which consider the marking 'correct'

	output
	------
	plot of 
	'''
	#import the data frames for the peak positions and the corresponding

	print ("This function may take a while as it compares all marked transits to the injected transits..")

	transit_locs_DF = transit_loc(data, WF)

	pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data_snr.txt', format='ascii',comment='#')

	for rown_n, row in data.iterrows(): # do for each row indivisually

		xmin = row['min_time']
		xmax = row['max_time']
		
		xminpix = row['min_pix']
		xmaxpix = row['max_pix']



		tic = row['candidate']
		
		pdf_peaks = row['pdf_peaks']  # the marked peaks. 
		transit_pix = row['x_peaks']
		
		transit_pix = string_to_list(transit_pix)
		pdf_peaks = string_to_list(pdf_peaks)


		transit_days = pixel_to_time(transit_pix,xmin, xmax, xminpix, xmaxpix)
		
		#transit_inj = transit_loc(tic) # the locations of the injected transits

		injtrans_pl = (transit_locs_DF.loc[tic,'injtrans_pl'])
		injtrans_eb = (transit_locs_DF.loc[tic,'injtrans_eb'])

		# for each TIC, loop through the injected transits and throught the marked transits and compare them
		# if any marked and injected are within 1 (?) day of one another, they are deemed as correctly marked
		
		if injtrans_pl != None:
			idx = []
			for i, t_mrk in enumerate(transit_days):
				for t_inj in injtrans_pl:
					if t_inj > (t_mrk - tolerance) and t_inj < (t_mrk + tolerance): # if the marked and injected planets are within one day of one another, plot it in green
						idx.append(i)

			for i in range(len(transit_days)):
				if i in idx:
					plt.scatter(transit_days[i], pdf_peaks[i], color = 'green')
				else:
					plt.scatter(transit_days[i], pdf_peaks[i], color = 'orange')
			
		if injtrans_eb != None:
			idx = []
			for i, t_mrk in enumerate(transit_days):
				for t_inj in injtrans_eb:
					if t_inj > t_mrk - tolerance and t_inj < t_mrk + tolerance:
						idx.append(i)
			
			for i in range(len(transit_days)):
				if i in idx:
					plt.scatter(transit_days[i], pdf_peaks[i], color = 'blue')
				else:
					plt.scatter(transit_days[i], pdf_peaks[i], color = 'orange')

		elif injtrans_pl == None and injtrans_eb == None:
			for i in range(len(transit_days)):
				plt.scatter(transit_days[i], pdf_peaks[i], color = 'red')
	
	plt.xlabel("Marked Peak (days)")
	plt.ylabel("Confidence")
	#plt.ylim(-0.3,2)
	plt.xlim(0,30)
	
	plt.savefig(outpath + '/marked_injected_comp_tol{}.png'.format(tolerance), format='png', dpi=300,bbox_inches='tight')






