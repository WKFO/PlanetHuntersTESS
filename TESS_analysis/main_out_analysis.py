#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import sys, os
import numpy as np
import pandas as pd 
import ujson

import matplotlib.pyplot as plt
from astropy.table import Table
from argparse import ArgumentParser
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as mpatches
import random

import gc
from argparse import ArgumentParser
from mpl_toolkits.mplot3d import Axes3D
from glob import glob

import out_analysis_cdpp
import ROC
import SNR
import class_count
import marking_location as ml


# e.g. python main_out_analysis.py 7865 52 --data-rls=ETE-6

active_workflow_id = 8349
active_workflow_major = 12  # subject to change

# what data release are we looking at e.g. ETE-6

outpath = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/figures{}'.format(active_workflow_id)

data_release = 'Rel01'

# import the data files
#classfile_in_WF = "/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/planet_aggregations_8349_12.csv"
#participants_WF = "/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/class_counts_8349_12.csv"

# THIS IMPORTS THE BETA TEST FILE FOR NOW- WILL CHANGE (FORMAT THE SAME SO GOOD FOR TESTING THE CODE)

classfile_in_WF = "/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/planet_aggregations_8312_9.csv"
participants_WF = "/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/class_counts_8349_12.csv"



# this is only for the ETE-6 data
pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data.txt', format='ascii',comment='#')

data_WF = pd.read_csv(classfile_in_WF)
who_WF = pd.read_csv(participants_WF)


#----------#----------#----------#----------#----------#----------#----------#

ml.marking_plot(data_WF, sort = 'p_Yes', dist = 'gauss', number = 3)  # doesn't plot anything at the moment.
ml.marking_method_comp(data_WF, sort = 'p_Yes', dist = 'gauss', number = 3)  # doesn't plot anything at the moment.

# dist options: all, box, gauss, DB
# sort options: #p_Yes, p_Yes_weight, pdf_peak_max,x_peak_max_box
ml.marking_plot(data_WF, sort = 'p_Yes', dist = 'box', number = 3)  # doesn't plot anything at the moment.
ml.marking_method_comp(data_WF, sort = 'p_Yes', dist = 'all', number = 3)  # doesn't plot anything at the moment.



# plots of the distributions of weights and number of classifications per user - just to get an idea of how well the project is going. 
def hist(data, variable):

	var = who_WF[variable]

	plt.hist(var,  bins=30)
	plt.xlabel("{}".format(variable))
	plt.ylabel("Frequency")
	
	plt.savefig(outpath + '/hist_{}_{}.png'.format(variable,active_workflow_id), format='png', dpi=80)

#hist(who_WF, "weight")
#hist(who_WF, "nclass_user")


# ---------------------------------------------------
#                      SIMS ONLY                    #
# ---------------------------------------------------

# Evaluate how well people are classifying the simulations

# therefore only feed in the simulated data

data_sim = data_WF[data_WF.subject_type  == True]

class_count.pie_identification_plot(data_sim)













