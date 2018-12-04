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


data_release = "Rel01"

active_workflow_ids = [8349]


active_workflow_major = [12]

colors = ['orange']


for i,active_workflow_id in enumerate(active_workflow_ids):
	participants_WF = "/Users/Nora/Documents/research/TESS/planethunters/{}/output/{}/class_counts_{}_{}.csv".format(data_release, active_workflow_id, active_workflow_id, active_workflow_major[i])

	who_WF = pd.read_csv(participants_WF)

	var = who_WF['nclass_user']
	
	#vals = np.histogram(a, bins=10, range=None, normed=False)

	plt.hist(var.dropna(),  bins=30, normed = True, alpha = 0.7, color=  colors[i], label = label[i])


plt.xlabel("Number of Classifications")
plt.ylabel("Frequency")

plt.legend()
plt.show()













