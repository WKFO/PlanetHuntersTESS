#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.table import Table
from glob import glob

ROOT = '/Users/Nora/Desktop/research/data/TESS/ETE-6/'
lcdir = ROOT + 'simulated_data/ALL'
pl = Table.read(ROOT + 'ground_truth/ete6_planet_data.txt', format='ascii',comment='#')
pl_tic = pl['col1'].flatten()
npl = len(pl_tic)
pl_depth = pl['col10'].flatten()
pl_dur = pl['col9'].flatten()
durs = np.array([0.5, 1, 2])
durs_ = np.array(['0_5', '1_0', '2_0'])
pl_cdpp = np.zeros(npl)
for i in range(npl):
    print(i,npl,pl_tic[i])
    lcfile = glob(lcdir + '/tess*{:d}*lc.fits'.format(pl_tic[i]))[0]
    hdul = pf.open(lcfile)
    hdr = hdul[1].header
    j = np.argmin(abs(pl_dur[i] - durs))
    pl_cdpp[i] = float(hdr['CDPP' + durs_[j]])
    hdul.close()
pl_snr = pl_depth / pl_cdpp
pl['col16'] = pl_snr
pl.write(ROOT + 'ground_truth/ete6_planet_data_snr_2.txt', format='ascii', overwrite = True)
