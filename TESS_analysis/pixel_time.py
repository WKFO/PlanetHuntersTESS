
import numpy as np
import pandas as pd 
from glob import glob
import astropy.io.fits as pf

from astropy.table import Table
import matplotlib.pyplot as plt
import ast

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


def pixel_to_time(xvals):
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
	time_per_pixel = (xmax - xmin) / (xmaxpix - xminpix)


	x_time = []

	for p in xvals:
		
		delta_pix = p - xminpix
	
		delta_time = delta_p * time_per_pixel
	
		x_t = xmin + delta_time

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











