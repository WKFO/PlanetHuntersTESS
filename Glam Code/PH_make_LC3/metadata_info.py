#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.table import Table
from glob import glob
import csv


ROOT = '/Users/Nora/Desktop/research/data/TESS/ETE-6/'
lcdir = ROOT + 'simulated_data/ALL'

pl = Table.read(ROOT + 'ground_truth/ete6_planet_data.txt', format='ascii',comment='#')
st = Table.read(ROOT + 'ground_truth/ete6_star_data.txt', format='ascii',comment='#')
eb = Table.read(ROOT + 'ground_truth/ete6_eb_data.txt', format='ascii',comment='#')

LC_type = [pl,st,eb]

tic_all = []
cdpp05_all = []
cdpp10_all = []
cdpp20_all = []

for LC in LC_type:
	tic = LC['col1'].flatten()

	n = len(tic)

	durs = np.array([0.5, 1, 2])
	durs_ = np.array(['0_5', '1_0', '2_0'])
	
	tics = np.zeros(n)
	cdpp05 = np.zeros(n)
	cdpp10 = np.zeros(n)
	cdpp20 = np.zeros(n)

	for i in range(n):
	    print(i,n,tic[i])
	    lcfile = glob(lcdir + '/tess*{:d}*lc.fits'.format(tic[i]))[0]
	    hdul = pf.open(lcfile)
	    hdr = hdul[1].header
	    cdpp05_all.append(float(hdr['CDPP0_5']))
	    cdpp10_all.append(float(hdr['CDPP1_0']))
	    cdpp20_all.append(float(hdr['CDPP2_0']))
	    tic_all.append(tic[i])
	    hdul.close()

# have a list of lists - mnake into one list 

with open(ROOT + 'ground_truth/ete6_all_metadata_info.txt', 'w') as f: # save in the photometry folder
    
    writer = csv.writer(f, delimiter=',')
    writer.writerows(zip(tic_all, cdpp05_all, cdpp10_all, cdpp20_all)) # zip the desired columns together


