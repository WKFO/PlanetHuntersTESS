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
import pixel_time_convert as ptc
import transit_loc
from itertools import chain, repeat, islice


directory = '/Users/Nora/Documents/research/TESS/Rel1'

# python main_simulated.py data simulated --binfac=7 --chunks=4 --no-binned --number=5 --min-snr=10

# This code makes 4(?) images out of each lightcurve and saves them as individual subjects. 
# This is so that the resolution is better for each image - easier to spot single transits (lots in TESS data)
# Previously there was one manifest input for multiple images - we now need one manifest input (line) for each.


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

    seed_num = 3 # the set of random samples #use seed =2 for the current images on zooniverse

    args = ap.parse_args()

    # import these tables - the planet table contains the SNR, which is good to make sure that the planets have an OK SNR
    pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data_snr.txt', format='ascii',comment='#')
    eb = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_eb_data.txt', format='ascii',comment='#')

    
    if args.min_snr == None: # THE CODE WILL FAIL...
        #lcfile0 = np.sort(glob('{}{:s}/tess*lc.fits'.format(directory,args.lcdir)))
        #nlc = len(lcfile0)
        #
        #if args.number != None: # randomly select the targets
        #    random.seed(seed_num)  # the ETE-6 data were all created with seed = 1. The images for examples with seed=7, 9
        #    random.shuffle(lcfile0)
        #    lcfile0 = lcfile0[:args.number]

        print ("to make simulated data you need to have a mininum SNR set...") 

    # get the tics of those planets (to inject) with SNR greater than a certain amount.
    else: 
        l = np.where(pl['col16'] > args.min_snr) # create a mask of elements where the SNR is greater than min
        tics_inj = np.unique(pl[l]['col1'].flatten())
        nlc = len(tics_inj)

        if args.number != None: # randomly select the planet targets to inject
            random.seed(seed_num)
            random.shuffle(tics_inj)
            tics_inj = tics_inj[:args.number]
            nlc = len(tics_inj)


    # also get random lightcurves from the real data - same number as tics
    lcfile0 = np.sort(glob('{}/{:s}/tess*lc.fits'.format(directory,args.lcdir)))  # will later want to filer out the ones where we know there are planets.
    random.seed(seed_num)  
    random.shuffle(lcfile0)  # randomly shuffle them (according to a seed - can change the seed but helps to have ome to keep track of what's going on.)
    lcfile0 = lcfile0[:len(tics_inj)]  # the same number of ranfom files as tics used to inject planets. 

    # -----------------------------
    # we now have a list of tic IDs - assuming that we entered a min SNR (otherwise we don't and this code will crash which is unfortunate)

    # find the chunks etc - Chunks needs to be the same as for the non sim data, otherwise that completely and utterly defeats the point of this excercise. 
    chunk = (args.tmax + args.overlap * (args.chunks-1)) / float(args.chunks)
    step = chunk - args.overlap
    sb.set(style="white", context="paper")


    # create a fun figure - again params need to be the same as for real data.
    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(12,5))
    plt.subplots_adjust(left=0.005, right=0.995, top=0.99, bottom=0.01)

    im = [ [] for i in range(args.chunks) ] # create the number of empty lists as there are chunks
    im_name  = []
    Magnitude = []
    Temperature = []
    Radius = []
    Target = []
    TICID = []  # main TIC ID is the ome of the planet - this will be used in the analysis code to pull out the informatino of where the transit are.
    Sector = []
    Web_link = []
    image = []
    TICLC = []  # store the tic ID Of the LC that was used to make the sims

    xmin_pix = []
    xmax_pix = []
    
    xmin_time = []
    xmax_time = []
    sim = []  # have a column to specify that this is a sim. Just to make sure that this is clear.


    # for the sim data we want some feedback - so make lots of empty lists to fill
    feedback_name = []
    feedback_id_name = []
    feedback_width_name = []
    feedback_success_name = []
    feedback_failure_name = []
    feedback_tolerance_name = []

    longest_list = 10 # for now assume that no LC will have more than 10 transits - may need to change
    
    feedback           = [ [] for i in range(longest_list) ]
    feedback_id        = [ [] for i in range(longest_list) ]
    feedback_width     = [ [] for i in range(longest_list) ]
    feedback_success   = [ [] for i in range(longest_list) ]
    feedback_failure   = [ [] for i in range(longest_list) ]
    feedback_tolerance = [ [] for i in range(longest_list) ]

    # start the real work. Loop through each tic id.
    for i in range(nlc):

        lcfile = lcfile0 # open one of the data lightcurves
        hdul   = pf.open(lcfile[i]) 

        # open the text file with the injected transit data. A file of mostly 0s and a couple of smaller numbers. 
        # first get the file name
        
        injected_pl = glob("/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/Injected/Planets/Planets_*{:d}.txt".format(tics_inj[i]))
        
        # now open it.
        inj_pl = np.genfromtxt(str(injected_pl[0]))        

        tic_LC     = int(hdul[0].header['TICID'])
        sec     = int(hdul[0].header['SECTOR'])
        cam     = int(hdul[0].header['CAMERA'])
        chi     = int(hdul[0].header['CCD'])
        tessmag = (hdul[0].header['TESSMAG']) # metadata 
        teff    = int(hdul[0].header['TEFF']) # metadata 
        srad    = (hdul[0].header['RADIUS'])  # metadata 

        scc     = '%02d%1d%1d' % (sec,cam,chi)
        
        print('TIC {:9d} SCC {:4s} ({:d} of {:d})'.format(tic_LC, scc, i+1, nlc))   # in the command line print how many files to go
        
        url = "https://archive.stsci.edu/tess/ete-6.html"  # will want to generate a link to the specific star but use this for now


        l_pl1 = np.where((pl['col1'] == tics_inj[i]) * (pl['col16'] >= args.min_snr))[0]        
        l_pl2 = np.where((pl['col1'] == tics_inj[i]) * (pl['col16'] < args.min_snr))[0]        
        l_eb = np.where(eb['col1'] == tics_inj[i])[0]

        data    = hdul[1].data
        time    = data['TIME']
        flux    = data['PDCSAP_FLUX']
        qual    = data['QUALITY']

        # --------------------------------------
        # before you do anything else with the flux of the LC - i.e. before binning and chunking etc, need to inject the transit.
        flux = inj_pl * flux
        
        # --------------------------------------

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

        # -------------------------------------------------------------------------
        # loop though each chunk, for each chunk create a line in the manifest and save one image
        # therefore for each TIC ID will have 4 manifest inputs
        # -------------------------------------------------------------------------

        for j in range(args.chunks):
            xmin = j * step
            xmax = xmin + chunk
            ls1 = (time1 >= xmin) * (time1 <= xmax)  # create a mask for the points that are not within xmin and xmax
            if not args.no_binned:
                ls2 = (time2 >= xmin) * (time2 <= xmax)
            ichunk = j + 1

            ax.cla()  # clear the figure before plotting more data onto it

            nm = 0
            transit_duration = []

            for l in l_pl1:
                P = pl['col3'][l]
                T0 = pl['col4'][l] - t0
                transit_duration = pl['col9'][l]
                while T0 < time1[ls1].min():
                    T0 += P
                while T0 < time1[ls1].max():
                    if not args.no_mark and args.min_snr != None:
                        plt.axvline(T0,color='chartreuse',lw=5,alpha=0.5)
                    nm += 1
                    T0 += P

            for l in l_pl2:
                P = pl['col3'][l]
                T0 = pl['col4'][l] - t0
                while T0 < time1[ls1].min():
                    T0 += P
                while T0 < time1[ls1].max():
                    if not args.no_mark and args.min_snr != None:
                        plt.axvline(T0,color='red',lw=5,alpha=0.5)
                    nm += 1
                    T0 += P

            for l in l_eb:
                P = eb['col3'][l]
                T0 = eb['col4'][l] - t0
                while T0 < time1[ls1].min():
                    T0 += P
                while T0 < time1[ls1].max():
                    if not args.no_mark and args.min_snr != None:
                        plt.axvline(T0,color='yellow',lw=5,alpha=0.5)
                    nm += 1
                    T0 += P

            if nm > 0:
                im_name.append('tess%09d_scc%04s_%02d_%09d_simulated.png' % (tic_LC,scc,ichunk,tics_inj[i])) # append the name of the image
    
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # PLOT the image
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
                points, = ax.plot(time1[ls1],flux1[ls1], 'w.', markersize=4) # needed to convert into pixels
    
                if not args.no_binned:
                    points, = ax.plot(time1[ls1],flux1[ls1],'co') # needed to convert into pixels
    
                # define that length on the x axis - don't want it to display the 0 point
                ax.set_xlim(xmin, xmax)  
    
                plt.subplots_adjust(left=0.001, right=0.998, top=0.999, bottom=0.005)
                # Get the x and y data and transform it into pixel coordinates
                x, y = points.get_data()
                xy_pixels = ax.transData.transform(np.vstack([x,y]).T)
                xpix, ypix = xy_pixels.T
    
                # In matplotlib, 0,0 is the lower left corner, whereas it's usually the upper 
                # right for most image software, so we'll flip the y-coords...
                xmin_time.append(xmin)
                xmax_time.append(xmax)
    
                xmin_pix.append(min(xpix))
                xmax_pix.append(max(xpix)) 
    
                delta_flux = np.max(flux1[ls1]) - np.min(flux1[ls1]) 
    
                percent_change = delta_flux * 0.1  
    
                ax.set_ylim(np.min(flux1[ls1]) - percent_change, np.max(flux1[ls1]) + percent_change)            # do we want to change this to be per figure and not per LC? 
                
                ax.xaxis.set_label_coords(0.063, 0.06) # position of the x-axis label 
    
    
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
    
                ax.set_axis_bgcolor("#03012d")
    
                # change how the image is saved.
    
                if args.chunks == 1:
                    print ("This code is not optimised for one chunk")
    
                else:
                    plt.savefig('%s/lightcurves/simulated/tess%09d_scc%04s_%02d_%09d_simulated.png' % (directory,tic_LC,scc,ichunk,tics_inj[i]), format='png')
        
                Magnitude.append("%.3f mag" % tessmag) 
                Temperature.append("{} K".format(teff))
                Radius.append("{} Solar Radii".format(str(srad)[0:5]))
                TICID.append(tics_inj[i])
                TICLC.append(tic_LC)
                Sector.append(sec)
                Web_link.append(url)
                sim.append(True)
    
                if tics_inj[i] in pl['col1']:
                    Target.append("planet")
                elif tics_inj[i] in eb['col1']:
                    Target.append("EB")
                else:
                    Target.append("star")
                
    
                injtrans_pl, injtrans_eb, tic = transit_loc.transit_loc_single(tics_inj[i], args.min_snr, xmin, xmax)
    
                transit_duration = ptc.deltat_to_pixel(transit_duration, xmin, xmax)

                ## find oversll longest list 
                #longest_list = 0
                #
                #for i in transit_locs_DF00['injtrans_pl']:
                #    if i != None:
                #        if len(i) > longest_list:
                #            longest_list = len(i)
    
                # select only the values that are in that chunk - for all chunks (have to loop through them again unofortunately)
    
                #transit_locs_DF0 = (xmin <= transit_locs_DF00['injtrans_pl']) & (transit_locs_DF00['injtrans_pl'] <= xmax)
                #transit_locs_DF = (xmin <= transit_locs_DF0['injtrans_eb']) & (transit_locs_DF0['injtrans_eb'] <= xmax)
    
    
                for j,trans_loc in enumerate(injtrans_pl):
                    
                    if trans_loc == None:
                        for k in range(0,longest_list):
                            feedback[k].append(None)
                            feedback_id[k].append(None)
                            feedback_width[k].append(None)    
                            feedback_success[k].append(None)  
                            feedback_failure[k].append(None)  
                            feedback_tolerance[k].append(None)
                    else:

                        trans_loc = [i - transit_duration/8 for i in trans_loc]
                        transit_list = transit_loc.pad(trans_loc, longest_list, '')
                
                        for l,t in enumerate(transit_list):
                            feedback[l].append(t)
    
                            if t != '':
                                feedback_id[l].append("1111{}".format(l))
                                feedback_width[l].append(transit_duration/4)    
                                feedback_success[l].append("Awesome, you correctly marked a simulated transit!")  
                                feedback_failure[l].append("Not quite, but please don't be discouraged, transits can be very difficult to spot. You can also check out the field guide for more guidance.")  
                                feedback_tolerance[l].append(transit_duration/4*0.4)
                            else:
                                feedback_id[l].append(None)
                                feedback_width[l].append(None)    
                                feedback_success[l].append(None)  
                                feedback_failure[l].append(None)  
                                feedback_tolerance[l].append(None)
            

            ## feedback_1_id #feedback_1_x   #feedback_1_width   #feedback_1_successMessage  #feedback_1_failureMessage  #feedback_1_tolerance
    for i in range (longest_list):
        feedback_name.append("#feedback_{}_x".format(i))
        feedback_id_name.append("#feedback_{}_id".format(i))
        feedback_width_name.append('#feedback_{}_width'.format(i))
        feedback_success_name.append('#feedback_{}_successMessage'.format(i))
        feedback_failure_name.append('#feedback_{}_failureMessage'.format(i))
        feedback_tolerance_name.append('#feedback_{}_tolerance'.format(i))


    if args.chunks == 1:
        print ("This code is not optimised for one chunk")
    
    else:

        meta_names = ["#Image", "Magnitude","Temperature", "Radius", "#Target", "!TIC ID", "Sector", "!Web_link", "#xmin_pix", "#xmax_pix", "#xmin_time", "#xmax_time", "#sim", '#TICLC']
        metadata_header = (meta_names + feedback_name + feedback_id_name + feedback_width_name + feedback_success_name + feedback_failure_name + feedback_tolerance_name)
        
        meta_data = [im_name, Magnitude,Temperature, Radius, Target, TICID, Sector, Web_link, xmin_pix, xmax_pix, xmin_time, xmax_time, sim, TICLC]
        metadata_data = meta_data + feedback + feedback_id + feedback_width + feedback_success + feedback_failure + feedback_tolerance


    with open('%s/lightcurves/simulated/manifest.csv' % (directory), 'w') as f: # save in the photometry folder
        writer = csv.writer(f, delimiter=',')
        writer.writerow(metadata_header)
        writer.writerows(zip(*metadata_data)) # zip the desired columns together
    
    plt.close()
        
        
        
        
        
        



