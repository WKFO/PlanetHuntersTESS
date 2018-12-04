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
    ap.add_argument('--min-snr', type=float, help='min. SNR of transit/eclipse to mark', default=None)
    ap.add_argument('--no-mark', action='store_true')
    ap.add_argument('--no-binned', action='store_true')
    ap.add_argument('--number', type=int, help='The number of randomly selected subjects (default all in outdir folder)', default=None)

    seed_num = 2 # the set of random samples # use seed =2 for the ones on zooniverse

    args = ap.parse_args()
    pl = Table.read('ETE-6/ground_truth/ete6_planet_data_snr.txt', format='ascii',comment='#')
    
    if args.min_snr == None:
        lcfile0 = np.sort(glob('{:s}/tess*lc.fits'.format(args.lcdir)))
        nlc = len(lcfile0)
            
        if args.number != None: # randomly select the planet targets
            random.seed(seed_num)  # the ETE-6 data were all created with seed = 1. The images for examples with seed=7, 9
            random.shuffle(lcfile0)
            lcfile0 = lcfile0[:args.number]
            nlc = len(lcfile0)

    else:
        l = np.where(pl['col16'] > args.min_snr) # create a mask of elements where the SNR is greater than min
        tics = np.unique(pl[l]['col1'].flatten())
        nlc = len(tics)


        if args.number != None: # randomly select the planet targets
            random.seed(seed_num)
            random.shuffle(tics)
            tics = tics[:args.number]
            nlc = len(tics)

    eb = Table.read('ETE-6/ground_truth/ete6_eb_data.txt', format='ascii',comment='#')
    beb = Table.read('ETE-6/ground_truth/ete6_backeb_data.txt', format='ascii',comment='#')
    
    chunk = (args.tmax + args.overlap * (args.chunks-1)) / float(args.chunks)
    step = chunk - args.overlap
    sb.set(style="white", context="paper")

    # create figure

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(12,5))


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


            if args.chunks == 1:
                im[j].append('tess%09d_scc%04s.png' % (tic,scc))
            else:
                im[j].append('tess%09d_scc%04s_chunk%02d.png' % (tic,scc,ichunk))



        # if the chunk is larger than one, we also need an image that is zoomed out and  has the divisions 
        # this image will need to have a binning factor twice as large as defined in the function 

        if args.chunks > 1:

            binf = args.binfac * 2 # binning factor for the ful frame is half of the indicidual zoomed in frames
            arg_chunks = 1

            chunk_1 = (args.tmax + args.overlap * (arg_chunks-1)) / float(arg_chunks)
            step_1 = chunk_1 - args.overlap

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
            teff    = int(hdul[0].header['TEFF']) # metadata 
            srad    = (hdul[0].header['RADIUS'])  # metadata 
    
            scc     = '%02d%1d%1d' % (sec,cam,chi)
            print('TIC {:9d} SCC {:4s} ({:d} of {:d}) full frame'.format(tic, scc, i+1, nlc))   # in the command line print how many files to go
            
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
    
            flux_list.append(flux1)

    for i in range (args.chunks):
        im_name.append("#Image{}".format(i))

    if args.chunks == 1:

        meta_names = ["Magnitude","Temperature", "Radius", "#Target", "TIC ID", "Sector", "Web_link", "#xmin_pix", "#xmax_pix"]
        metadata_header = (im_name + meta_names)
    
        meta_data = [Magnitude,Temperature, Radius, Target, TICID, Sector, Web_link, xmin_pix, xmax_pix]
        metadata_data = im + meta_data
    
    else:
        meta_names = ["#Image", "Magnitude","Temperature", "Radius", "#Target", "TIC ID", "Sector", "Web_link", "#xmin_pix", "#xmax_pix"]
        metadata_header = (im_name + meta_names)
    
        meta_data = [image, Magnitude,Temperature, Radius, Target, TICID, Sector, Web_link, xmin_pix, xmax_pix]
        metadata_data = im + meta_data
    

    with open('ETE-6/Chunks=%s/%s/ML_tess.csv' % (args.chunks,args.outdir), 'w') as f: # save in the photometry folder
        writer = csv.writer(f, delimiter=',')
        writer.writerow(metadata_header)
        writer.writerows(zip(flux1, )) # zip the desired columns together
    
    plt.close()









