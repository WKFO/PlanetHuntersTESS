#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from astropy.table import Table
import seaborn as sb
from glob import glob
from argparse import ArgumentParser
from matplotlib.ticker import AutoMinorLocator
import csv
import matplotlib.patches as mpatches
import random
from os.path import basename

import pixel_time_convert as ptc
import transit_loc
from itertools import chain, repeat, islice

from time import time as clock

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    mpi_rank = comm.Get_rank()
    mpi_size = comm.Get_size()
    with_mpi = True
except ImportError:
    mpi_rank = 0
    mpi_size = 1
    with_mpi = False
 
mpi_root = 0
# np.seterr(all='ignore')

#in dir/Users/Nora/Documents/research/TESS/Rel1/data
# outdir = /Users/Nora/Documents/research/TESS/Rel1/lightcurves/glam_test


# python main_glam.py --lcdir=dir/Users/Nora/Documents/research/TESS/Rel1/data outdir=/Users/Nora/Documents/research/TESS/Rel1/lightcurves/glam_test --binfac=7 --chunks=4 --no-binned


#####################################################
# This is a standin until we get astroquery working #
#####################################################
url = "https://archive.stsci.edu/tess/ete-6.html"  # will want to generate a link to the specific star but use this for now

# python main.py data data --binfac=7 --chunks=4 --no-binned 

# This code makes 4(?) images out of each lightcurve and saves them as individual subjects. 
# This is so that the resolution is better for each image - easier to spot single transits (lots in TESS data)
# Previously there was one manifest input for multiple images - we now need one manifest input (line) for each.

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
                 new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def process_lc(lcfile, args):

    chunk = (args.tmax + args.overlap * (args.chunks-1)) / float(args.chunks)
    step = chunk - args.overlap

    # create figure
    sb.set(style="white", context="paper")
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(12,5))
    plt.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.01)



    hdul    = pf.open(lcfile)

    tic     = int(hdul[0].header['TICID'])
    sec     = int(hdul[0].header['SECTOR'])
    cam     = int(hdul[0].header['CAMERA'])
    chi     = int(hdul[0].header['CCD'])
    tessmag = hdul[0].header['TESSMAG']
    teff    = hdul[0].header['TEFF']
    srad    = hdul[0].header['RADIUS']

    scc     = '%02d%1d%1d' % (sec,cam,chi)
    
    print('TIC {:9d} SCC {:4s}'.format(tic, scc))   # in the command line print how many files to go
    
    data    = hdul[1].data
    time    = data['TIME']
    flux    = data['PDCSAP_FLUX']
    qual    = data['QUALITY']
    hdul.close()
    t0      = time[0]
    time    -= t0
    l       = np.isfinite(time) * np.isfinite(flux) * (qual == 0)
    time    = time[l]
    flux    = flux[l]
    m       = np.median(flux)
    flux    = 1e2 * ((flux / m) - 1)
    N       = len(time)

    n       = int(np.floor(N/args.binfac)*args.binfac)
    X       = np.zeros((2,n))
    X[0,:]  = time[:n]
    X[1,:]  = flux[:n]
    Xb      = rebin(X, (2,int(n/args.binfac)))
    time1    = Xb[0]
    flux1    = Xb[1]

    if not args.no_binned:
        n = int(np.floor(N/args.binfac/10)*args.binfac*10)
        X = np.zeros((2,n))
        X[0,:] = time[:n]
        X[1,:] = flux[:n]
        Xb = rebin(X, (2,int(n/args.binfac/10)))
        time2 = Xb[0]
        flux2 = Xb[1]


    # -------------------------------------------------------------------------
    # loop though each chunk, for each chunk create a new row in the manifest and save one image
    # i.e. each TIC ID will have 4 manifest inputs
    # -------------------------------------------------------------------------
    
    im_name = []
    xmin_pix = []
    xmax_pix = []    
    xmin_time = []
    xmax_time = []

    for j in range(args.chunks):
        xmin = j * step
        xmax = xmin + chunk
        ls1 = (time1 >= xmin) * (time1 <= xmax)  # create a mask for the points that are not within xmin and xmax
        if not args.no_binned:
            ls2 = (time2 >= xmin) * (time2 <= xmax)
        ichunk = j + 1

        ax.cla()  # clear the figure before plotting more data onto it

        # save the file name
        im_name.append('tess%09d_scc%04s_%02d.png' % (tic,scc,ichunk)) # append the name of the image


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # PLOT the image
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        points, = ax.plot(time1[ls1],flux1[ls1], 'w.', markersize=4) # needed to convert into pixels

        if not args.no_binned:
            points, = ax.plot(time1[ls1],flux1[ls1],'co') # needed to convert into pixels

        plt.subplots_adjust(left=0.001, right=0.998, top=0.999, bottom=0.005)

        # define that length on the x axis - don't want it to display the 0 point
        ax.set_xlim(xmin, xmax)  

        # Get the x and y data and transform it into pixel coordinates for the manifest (not sure how useful these will be but store anyway just in case_)
        x, y = points.get_data()
        xy_pixels = ax.transData.transform(np.vstack([x,y]).T)
        xpix, ypix = xy_pixels.T

        xmin_time.append(xmin)
        xmax_time.append(xmax)

        xmin_pix.append(min(xpix))
        xmax_pix.append(max(xpix)) 

        delta_flux = np.max(flux1[ls1]) - np.min(flux1[ls1]) 

        # set the y lim. 
        percent_change = delta_flux * 0.1  
        ax.set_ylim(np.min(flux1[ls1]) - percent_change, np.max(flux1[ls1]) + percent_change)  # do we want to change this to be per figure and not per LC? 
        
        # label the axis.
        ax.xaxis.set_label_coords(0.063, 0.06) # position of the x-axis label 


        # tick mark/axis parameters

        minorLocator = AutoMinorLocator()
        ax.xaxis.set_minor_locator(minorLocator)
        ax.tick_params(direction='in', which ='minor', colors='grey',length=3, labelsize=13)

        minorLocator = AutoMinorLocator()
        ax.yaxis.set_minor_locator(minorLocator)
        ax.tick_params(direction='in', length = 3, which ='minor', colors='grey', labelsize=13)

        ax.tick_params(axis="y",direction="in", pad= -20)
        ax.tick_params(axis="x",direction="in", pad= -17)   
        ax.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on')

        ax.set_xlabel("DAYS",fontsize = 10)
        #ax.artists.append(circle1)
        #ax.artists.append(circle2)

        # ax.set_axis_bgcolor("#03012d")
        ax.set_facecolor("#03012d")

        # save the image - make sure the name of the saved image is the same as the name in the manifest.

        plt.savefig('%s/%s' % (args.outdir,im_name[-1]), format='png')


    mnd = {}
    mnd['Magnitude'] = '{:.2f}'.format(float(tessmag))
    try:
        mnd['Temperature'] = '{:d}'.format(int(teff))
    except:
        mnd['Temperature'] = None
    try:
        mnd['Radius'] = '{:.3f}'.format(float(srad))
    except:
        mnd['Radius'] = None
    mnd['TICID'] = tic
    mnd['Sector'] = sec
    mnd['Web_link'] = url
    mnd['sim'] = False
    mnd['Target'] = False
    mnd['im_name'] = im_name
    mnd['xmin_time'] = xmin_time
    mnd['xmax_time'] = xmax_time
    mnd['xmin_pix'] = xmin_pix
    mnd['xmax_pix'] = xmax_pix

    return mnd        

if __name__ == '__main__':
    ap = ArgumentParser(description='Script to pLot TESS LCs for Zooniverse project')
    ap.add_argument('--lcdir', type=str, help='Directory containing input light curve files.', 
        default='/kepler/kepler2/TESS/ETE-6/light_curves/two_min/')
    ap.add_argument('--outdir', type=str, help='Directory to place output PNG files in.',
        default='/kepler/kepler2/TESS/planethunters/Rel00/')
    ap.add_argument('--binfac', type=int, help='Binning factor', default=3)
    ap.add_argument('--chunks', type=int, help='Chunks per sector', default=1)
    ap.add_argument('--overlap', type=float, help='Overlap between chunks (in days)', default=1.0)
    ap.add_argument('--tmax', type=float, help='Total duration of sector', default=27.8)
    ap.add_argument('--no-binned', action='store_true')
    ap.add_argument('--number', type=int, help='The number of randomly selected subjects (default all in outdir folder)', default=None)

    seed_num = 2 # the set of random samples #use seed =2 for the current images on zooniverse

    args = ap.parse_args()


    ## Master node
    ## -----------
    if mpi_rank == mpi_root:
        t0 = clock()

        # select the lightcurves
        lcfile0 = sorted(glob('{:s}/tess*lc.fits'.format(args.lcdir)))
        nlc = len(lcfile0)
    
        if args.number != None: # randomly select the targets
            random.seed(seed_num)  # the ETE-6 data were all created with seed = 1. The images for examples with seed=7, 9
            random.shuffle(lcfile0)
            lcfile0 = lcfile0[:args.number]

        metadata_header = ["#Image", "Magnitude","Temperature (K)", "Radius (solar radii)", "#Target", "!TIC ID", "Sector", "!Web_link", "#xmin_pix", "#xmax_pix", "#xmin_time", "#xmax_time", "#sim"]
    
        with open('%s/manifest.csv' % (args.outdir), 'w') as f: # save in the photometry folder
            writer = csv.writer(f, delimiter=',')
            writer.writerow(metadata_header)

        ## Without MPI or running with a single node
        ## =========================================
        if (not with_mpi) or (mpi_size==1) or (nlc==1):
            print ('processing {:d} LCs'.format(nlc))
            for f in lcfile0:
                print ('processing TIC {}'.format(int(basename(f).split('-')[1])))
                res = process_lc(f, args)
                with open('%s/manifest.csv' % (args.outdir), 'a') as f: # save in the photometry folder
                    writer = csv.writer(f, delimiter=',')
                    for ichunk in range(args.chunks):
                        metadata_data = [res['im_name'][ichunk]]
                        metadata_data.append(res['Magnitude'])
                        metadata_data.append(res['Temperature'])
                        metadata_data.append(res['Radius'])
                        metadata_data.append(res['Target'])
                        metadata_data.append(res['TICID'])
                        metadata_data.append(res['Sector'])
                        metadata_data.append(res['Web_link'])
                        metadata_data.append(res['xmin_pix'][ichunk])
                        metadata_data.append(res['xmax_pix'][ichunk])
                        metadata_data.append(res['xmin_time'][ichunk])
                        metadata_data.append(res['xmax_time'][ichunk])
                        metadata_data.append(res['sim'])
                        writer.writerow(metadata_data) 

        else:
            ## Master node
            ## -----------
            if mpi_rank == 0:
                free_workers = list(range(1,mpi_size))
                active_workers = []
                n_finished_items = 0

                while lcfile0 or active_workers:
                    ## Send a file
                    while lcfile0 and free_workers:
                        w = free_workers.pop()
                        f = lcfile0.pop()
                        print ('Sending TIC {} to worker {}'.format(int(basename(f).split('-')[1]), w))
                        comm.send(f, dest=w, tag=0)
                        active_workers.append(w)
                        
                        
                    ## Receive the results
                    for w in active_workers:
                        if comm.Iprobe(w, 2):
                            res = comm.recv(source=w, tag=2)

                            print ('Worker {} finished processing TIC {}'.format(w, res['TICID']))
                            with open('%s/manifest.csv' % (args.outdir), 'a') as f: # save in the photometry folder
                                writer = csv.writer(f, delimiter=',')
                                for ichunk in range(args.chunks):
                                    metadata_data = [res['im_name'][ichunk]]
                                    metadata_data.append(res['Magnitude'])
                                    metadata_data.append(res['Temperature'])
                                    metadata_data.append(res['Radius'])
                                    metadata_data.append(res['Target'])
                                    metadata_data.append(res['TICID'])
                                    metadata_data.append(res['Sector'])
                                    metadata_data.append(res['Web_link'])
                                    metadata_data.append(res['xmin_pix'][ichunk])
                                    metadata_data.append(res['xmax_pix'][ichunk])
                                    metadata_data.append(res['xmin_time'][ichunk])
                                    metadata_data.append(res['xmax_time'][ichunk])
                                    metadata_data.append(res['sim'])
                                    writer.writerow(metadata_data) 
                            free_workers.append(w)
                            active_workers.remove(w)
                            n_finished_items += 1

                print ('Completed {} of {} files'.format(n_finished_items,nlc))
                if n_finished_items < nlc:
                    print ('Failed on {} files'.format(nlc-n_finished_items))
                for w in free_workers:
                    comm.send(-1, dest=w, tag=0)
                t1 = clock()
                print('Time taken: {}'.format(t1-t0))

    ## Worker node
    ## -----------
    else:
        while True:
            filename = comm.recv(source=mpi_root, tag=0)
            if filename == -1:
                break
            res = process_lc(filename, args)
            comm.send(res, dest=mpi_root, tag=2)    

