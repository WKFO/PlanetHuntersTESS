#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import numpy as np
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
import matplotlib.patches as patches
import pandas as pd

import pixel_time_convert as ptc

# ETE-6/data/LC_ALL ETE-6/img_test --binfac=14 --chunks=2 --no-binned --number=1
directory = '/Users/Nora/Documents/research/TESS/'

def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
                 new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

if __name__ == '__main__':
    ap = ArgumentParser(description='Script to pLot TESS LCs for Zooniverse project')
    ap.add_argument('lcdir', type=str, help='Directory containing input light curve files.')
    ap.add_argument('outdir', type=str, help='Directory to place output PNG files in.')
    ap.add_argument('--binfac', type=int, help='Binning factor', default=3)
    ap.add_argument('--chunks', type=int, help='Chunks per sector', default=1)
    ap.add_argument('--overlap', type=float, help='Overlap between chunks (in days)', default=1.0)
    ap.add_argument('--tmax', type=float, help='Total duration of sector', default=27.8)
    ap.add_argument('--min-snr', type=float, help='min. threshold of marked transit', default=None)
    ap.add_argument('--no-mark', action='store_true')
    ap.add_argument('--no-binned', action='store_true')
    ap.add_argument('--number', type=int, help='The number of randomly selected subjects (default all in outdir folder)', default=None)

    seed_num = 2 # the set of random samples # use seed =2 for the ones on zooniverse

    args = ap.parse_args()
    
    pl = Table.read('{}ETE-6/ground_truth/ete6_planet_data_snr.txt'.format(directory), format='ascii',comment='#')
    
    # the table contains the peak positions of the gaussians marked in the initial version of the project.
    # for this puprpose we care about pdf_peaks and x_peaks

    df = pd.read_csv('{}ETE-6/output/WF1/planet_aggregations_7865_52.csv'.format(directory)) 
    
    for i, row in df.iterrows(): # do for each row indivisually

        pdf_peaks = row['pdf_peaks']
        x_peaks = row['x_peaks']
        
        x_peaks = ptc.string_to_list(x_peaks)
        pdf_peaks = ptc.string_to_list(pdf_peaks)
        x_hours = ptc.pixel_to_time(x_peaks)

        df.at[i, 'x_peaks'] = x_hours

    #pdf_peaks = string_to_list(pdf_peaks)
    #transit_days = pixel_to_time(transit_pix)

    TIC_wanted = df['TIC_ID']
    
    if args.min_snr == None:
        lcfile0 = []
        for TIC_ID in TIC_wanted:
            file = glob('{}{:s}/tess*{}*lc.fits'.format(directory,args.lcdir,TIC_ID))
            lcfile0.append(file)
            #print (file)
            #nlc = len(lcfile0)
    
            if args.number != None: # randomly select the planet targets
                random.seed(seed_num)  # the ETE-6 data were all created with seed = 1. The images for examples with seed=7, 9
                random.shuffle(lcfile0)
                lcfile0 = lcfile0[:args.number]

        nlc = len(lcfile0)

    else:
        l = np.where(clipping['pdf_peaks'] > args.min_snr) # create a mask of elements where the SNR is greater than min
        tics = np.unique(clipping[l]['pdf_peaks'].flatten())
        nlc = len(tics)

        if args.number != None: # randomly select the planet targets
            random.seed(seed_num)
            random.shuffle(tics)
            tics = tics[:args.number]
            nlc = len(tics)

    lcfile0 = [val for sublist in lcfile0 for val in sublist]

    eb = Table.read('{}ETE-6/ground_truth/ete6_eb_data.txt'.format(directory), format='ascii',comment='#')
    beb = Table.read('{}ETE-6/ground_truth/ete6_backeb_data.txt'.format(directory), format='ascii',comment='#')
    
    chunk = (args.tmax + args.overlap * (args.chunks-1)) / float(args.chunks)
    step = chunk - args.overlap
    sb.set(style="white", context="paper")

    # create figure
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8,8))

    im = [ [] for i in range(args.chunks) ] # create the number of empty lists as there are chunks
    im_name  = []
    Magnitude = []
    Temperature = []
    Radius = []
    Target = ["{}".format(args.lcdir[9:])] * nlc
    TICID = []
    Sector = []
    Web_link = []
    image = []
    xmin_pix = []
    xmax_pix = []
    im_num = []

    for i in range(nlc):

        if args.min_snr == None:
            lcfile = lcfile0
            hdul    = pf.open(lcfile[i])
        else:
            lcfile  = glob(args.lcdir + '/tess*{:d}*.fits'.format(tics[i]))

            if len(lcfile)==0:
                continue
            lcfile = lcfile[0]
            hdul    = pf.open(lcfile)

        tic     = int(hdul[0].header['TICID'])
        sec     = int(hdul[0].header['SECTOR'])
        cam     = int(hdul[0].header['CAMERA'])
        chi     = int(hdul[0].header['CCD'])
        tessmag = (hdul[0].header['TESSMAG']) # metadata 
        teff    = int(hdul[0].header['TEFF']) # metadata 
        srad    = (hdul[0].header['RADIUS'])  # metadata 

        scc     = '%02d%1d%1d' % (sec,cam,chi)
        
        print('TIC {:9d} SCC {:4s} ({:d} of {:d})'.format(tic, scc, i+1, nlc))   # in the command line print how many files to go
        
        url = "https://archive.stsci.edu/tess/ete-6.html"  # will want to generate a link to the specific star but use this for now

        if not args.no_mark and args.min_snr != None:
            l_pl1 = np.where((pl['col1'] == tic) * (pl['col16'] >= args.min_snr))[0]        
            l_pl2 = np.where((pl['col1'] == tic) * (pl['col16'] < args.min_snr))[0]        
            l_eb = np.where(eb['col1'] == tic)[0]
            l_beb = np.where(beb['col1'] == tic)[0]

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
        flux    = 1e3 * ((flux / m) - 1)
        N       = len(time)

        n       = int(np.floor(N/args.binfac)*args.binfac)
        X       = np.zeros((2,n))
        X[0,:]  = time[:n]
        X[1,:]  = flux[:n]
        Xb      = rebin(X, (2,int(n/args.binfac)))
        time1   = Xb[0]
        flux1   = Xb[1]

        if not args.no_binned:
            n = int(np.floor(N/args.binfac/10)*args.binfac*10)
            X = np.zeros((2,n))
            X[0,:] = time[:n]
            X[1,:] = flux[:n]
            Xb = rebin(X, (2,int(n/args.binfac/10)))
            time2 = Xb[0]
            flux2 = Xb[1]

        # create empty lists for the manifest file 

        Magnitude.append("%.3f mag" % tessmag) 
        Temperature.append("{} K".format(teff))
        Radius.append("{} Solar Radii".format(str(srad)[0:5]))
        TICID.append(tic)
        Sector.append(sec)
        Web_link.append(url)

        for j in range(args.chunks):
            xmin = j * step
            xmax = xmin + chunk
            ls1 = (time1 >= xmin) * (time1 <= xmax)
            if not args.no_binned:
                ls2 = (time2 >= xmin) * (time2 <= xmax)
            ichunk = j + 1

            ax.cla()  # clear the figure before plotting more data onto it

            if not args.no_mark and args.min_snr != None:
                nm = 0
                for l in l_pl1:
                    P = pl['col3'][l]
                    T0 = pl['col4'][l] - t0
                    while T0 < time1[ls1].min():
                        T0 += P
                    while T0 < time1[ls1].max():
                        plt.axvline(T0,color='chartreuse',lw=5,alpha=0.5)
                        nm += 1
                        T0 += P
                for l in l_pl2:
                    P = pl['col3'][l]
                    T0 = pl['col4'][l] - t0
                    while T0 < time1[ls1].min():
                        T0 += P
                    while T0 < time1[ls1].max():
                        plt.axvline(T0,color='red',lw=5,alpha=0.5)
                        nm += 1
                        T0 += P
                for l in l_eb:
                    P = eb['col3'][l]
                    T0 = eb['col4'][l] - t0
                    while T0 < time1[ls1].min():
                        T0 += P
                    while T0 < time1[ls1].max():
                        plt.axvline(T0,color='yellow',lw=5,alpha=0.5)
                        nm += 1
                        T0 += P
                for l in l_beb:
                    P = beb['col3'][l]
                    T0 = beb['col4'][l] - t0
                    while T0 < time1[ls1].min():
                        T0 += P
                    while T0 < time1[ls1].max():
                        nm += 1
                        plt.axvline(T0,color='orange',lw=5,alpha=0.5)
                        T0 += P


            # plotting

            tic_df = df.loc[df['TIC_ID'] == tic]
            peak_list = list(tic_df['x_peaks'])[0]

            for i,peak in enumerate(peak_list):

                im[j].append('tess%09d_scc%04s_%s.png' % (tic,scc,i))

                ax.plot(time1[ls1],flux1[ls1],'w.', markersize=4)
    
                if not args.no_binned:
                    ax.plot(time2[ls2],flux2[ls2],'co')
    
                # define that length go the x axis - don't want it to display the 0 point and want it to show 30 days

                ax.set_xlim(peak - 1, peak + 1)

                ax.set_ylim(np.min(flux1) - 0.5, np.max(flux1) + 0.5)
                ax.xaxis.set_label_coords(0.085, 0.06) # position of the x-axis label 
    
                minorLocator = AutoMinorLocator()
                ax.xaxis.set_minor_locator(minorLocator)
                ax.tick_params(direction='in', which ='minor', colors='w',length=3, labelsize = 13)
                
                minorLocator = AutoMinorLocator()
                ax.yaxis.set_minor_locator(minorLocator)
                ax.tick_params(direction='in', length = 3, which ='minor', colors='w', labelsize = 13)
                
                ax.tick_params(axis="y",direction="in", pad= -20)
                ax.tick_params(axis="x",direction="in", pad= -17)   
                ax.tick_params(axis='both', length = 7, left='on', top='on', right='on', bottom='on', labelsize = 13)
    
                ax.set_xlabel("DAYS", fontsize = 16)
                #ax.artists.append(circle1)
                #ax.artists.append(circle2)
                
                ax.set_axis_bgcolor("#03012d")
                
                im_num.append(i)
                # change how the image is saved.
                
                plt.savefig('%sETE-6/Chunks=%s/%s/tess%09d_scc%04s_%s.png' % (directory,args.chunks,args.outdir,tic,scc,i), format='png',dpi = 80, bbox_inches='tight')
    

                plt.tight_layout()
                xmin_pix.append((ax.transData.transform((0, 0))[0]-2)*0.798319) # this number is a scaling factor needed for the simulated data (change the code for real data to not need this)
                xmax_pix.append((ax.transData.transform((xmax, 0))[0]-6)*0.798319) 
    

    if args.chunks == 1:

        meta_names = ["Image", "Magnitude","Temperature", "Radius", "#Target", "!TIC ID", "Sector", "!Web_link", "#xmin_pix", "#xmax_pix", "#im_num"]
        metadata_header = (meta_names)
    
        meta_data = [Magnitude,Temperature, Radius, Target, TICID, Sector, Web_link, xmin_pix, xmax_pix, im_num]
        metadata_data = im + meta_data
    
    
    with open('%sETE-6/Chunks=%s/%s/manifest.csv' % (directory,args.chunks,args.outdir), 'w') as f: # save in the photometry folder
        writer = csv.writer(f, delimiter=',')
        writer.writerow(metadata_header)
        writer.writerows(zip(*metadata_data)) # zip the desired columns together
    
    plt.close()




