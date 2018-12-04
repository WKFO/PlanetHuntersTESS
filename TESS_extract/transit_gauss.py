import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import get_info
import peakutils
from astropy.table import Table
from glob import glob
import astropy.io.fits as pf


def xval_hist_box(grp):
    
    """
    
    Looks at all of the answers for one subject (i.e. all the markings for one lightcurve)
    Returns the net (summed and weighted) distribution of answers

    Args:
        grp_tic (pd.DataFrame): for one TIC only (pre-filtered), classifications by many users

    Returns:
        (list): points at which distribution is evaluated
        (list): distribution of answers, summed and weighted
    """

    grp = grp.dropna(subset=['xvals', 'width']) # get rid of nan, None etc values (these will disrupt the code)

    min_pixels = 0
    max_pixels = 1193

    x_eval = np.linspace(min_pixels, max_pixels, max_pixels)
    grp_pdf = np.zeros_like(x_eval)

    for row_n, row in grp.iterrows():
        xvals = row['xvals']
        widths = row['width']
        weight = row['weight']

        try:
            #user_vote_pdf = get_info.get_net_voting_distribution(xvals, widths, weight, x_eval)
            #grp_pdf += user_vote_pdf
            
            user_vote_pdf = get_info.get_net_voting_distribution_box(xvals, widths, weight, x_eval)
            grp_pdf += user_vote_pdf

        except:
            print ("no marking")

    return x_eval, grp_pdf


def xval_hist_gauss(grp):
    
    """
    
    Looks at all of the answers for one subject (i.e. all the markings for one lightcurve)
    Returns the net (summed and weighted) distribution of answers

    Args:
        grp_tic (pd.DataFrame): for one TIC only (pre-filtered), classifications by many users

    Returns:
        (list): points at which distribution is evaluated
        (list): distribution of answers, summed and weighted
    """

    grp = grp.dropna(subset=['xvals', 'width']) # get rid of nan, None etc values (these will disrupt the code)

    min_pixels = 0
    max_pixels = 1180

    x_eval = np.linspace(min_pixels, max_pixels, 1180)
    grp_pdf = np.zeros_like(x_eval)

    for row_n, row in grp.iterrows():
        xvals = row['xvals']
        widths = row['width']
        weight = row['weight']

        try:
            user_vote_pdf = get_info.get_net_voting_distribution(xvals, widths, weight, x_eval)
            grp_pdf += user_vote_pdf

        except:
            print ("no marking")

    return x_eval, grp_pdf



def get_distribution_peaks(x_eval, grp_pdf):
    """
    Function to return the x values and tic pdf of the peaks. 

    # record the value of the pdf peaks
    # record the (float) pixels with a pdf peak

    # TODO check for what happens if NO-ONE answers
    """

    # find indices of distribution peak
    indexes = peakutils.indexes(grp_pdf, thres=0.05, min_dist=2)

    if len(indexes) == 0:
        return np.zeros(1),np.zeros(1)
   
    else:

        return x_eval[indexes], grp_pdf[indexes]


def get_max_peak(x_eval, grp_pdf):
    """
    Function to return the x values and tic pdf of the highest peak
    """

    # find indices of the highest point - i.e. where grp_pdf is max

    indexes = peakutils.indexes(grp_pdf, thres=0.05, min_dist=2)
    if len(indexes) == 0:
        return 0,0
   
    else:
        peaks_x = list(x_eval[indexes])
        peaks_pdf = list(grp_pdf[indexes])

        if len(peaks_x) == 1:
            return float(peaks_x[0]), float(peaks_pdf[0])
        
        elif len(peaks_x) > 1:
            return float(max(list(peaks_x))), float(max(list(peaks_pdf)))


    # # plot everything
    # fig, (ax1, ax2) = plt.subplots(2,1, figsize =(10,10))  # define the axes to plot
    # #get a gaussian distribution of probabilities for each marked transit
    # im = mpimg.imread("/Users/Nora/Documents/research/TESS/ETE-6/Chunks=1/zooniverse_upload_100/{}".format(grp_tic['candidate'].unique()[0])) 
    # ax1.imshow(im)
    
    # #x_eval = [(x - 45.89) * 0.039 for x in x_eval]
    # ax2.plot(x_eval, tic_pdf)
    # ax2.plot(x_eval[indexes], tic_pdf[indexes], marker= 'x',color = 'r', ls='none')
    # ax2.set_xlabel('Pixels')
    # ax2.set_ylabel('Net PDF from all users, weighted')
    # plt.savefig('/Users/Nora/Documents/research/TESS/ETE-6/output/results_figs/output_PDM/hist{}'.format(grp_tic['candidate'].unique()[0]), format='png', bbox_inches='tight')
    
    # ax1.cla()
    # ax2.cla()
    

pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data.txt', format='ascii',comment='#')
eb = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_eb_data.txt', format='ascii',comment='#')

# combined list of TICS with a transit



def transit_loc_per_tic(tic, min_time, max_time):

    '''

    A function to find the times of the injected transits for a given TIC ID (both planets and EB). requires the ground truth file. 
    This will not be relevant for the real project as we will have simulated data in there. Only necessary on simulated data. 

    input:
    -----
    tic (int)

    output:
    ------
    transit_loc_TIC (list): list of all the transit times for that TIC ID. If no transits are found i.e. there is no planet, not EB, return [0,0] - these can later be ignored.

    '''

    # for this TIC, where are the transits?

    # combine the lists of the EB and the planets - we count both as a transit (might change later?)
    tics_pl_eb  = list(pl['col1']) + list(eb['col1'])  #tic
    P_pl_eb = list(pl['col3']) + list(eb['col3'])  # period
    T0_pl_eb = list(pl['col4']) + list(eb['col4'])  # initial time
    # eb's do not have a duration... great.

    width_pl_eb0 = list(pl['col9']) + ([1.5] * len(eb['col4']))   # eb's do not have a duration... great.

    width_pl_eb = [i/24 for i in width_pl_eb0]

    if float(tic) in tics_pl_eb:  #if subject has a planet or EB then...

        # find the fits file for that tic (this method is slow, I will want to create a seperate file with all this information in so that it doens't have to run everytime)
        lcfile  = glob("/Users/Nora/Documents/research/TESS/ETE-6/data/LC_ALL/tess*{}*.fits".format(tic))
        
        # open the file and get the header information
        lcfile = lcfile[0]
        hdul    = pf.open(lcfile)
    
        data    = hdul[1].data
        flux    = data['PDCSAP_FLUX']
        qual    = data['QUALITY']
        time    = data['TIME']
        hdul.close()

        t0      = time[0]
        time    -= t0
        l       = np.isfinite(time) * np.isfinite(flux) * (qual == 0)
        time    = time[l]
        
        l = tics_pl_eb.index(float(tic))
        #l_pl_eb = np.where((tics_pl_eb == tic))[0] 

        transit_loc_TIC = []
        widths = []

        P = P_pl_eb[l]
        widths.append(width_pl_eb[l])  # get the width
        T0 = T0_pl_eb[l] - t0

        while T0 < min_time:
            T0 += P
        while T0 < max_time:
            transit_loc_TIC.append(T0)
            T0 += P
        
        # a list of the transits for that tic and a list of a single value of the width of all those transits.
        return transit_loc_TIC, widths 

    else:
        return [0,0], [0]



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



def deltat_to_pixel(deltatime, tmin, tmax):
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
    delta_pix = []

    for t in deltatime:
        x_pixel = (t / (tmax - tmin)) * 1193 # not sure why this is 
        delta_pix.append(x_pixel)

    return delta_pix



def peaks_right_wrong(grp):
    '''
    input is all the classifications fot one tic ID i.e. a pandas array of the data for one TIC ID 
    this means that I can iimport the ground truth for that one TIC_ID

    I then look at each row (classificatiom) individually and assess whether the marked box overlaps with the marked transit
    
    A correctly drawn box over the transit gets a point of 1, otherwise 0.

    For each transit for each TIC return the number of correctly drawn boxes - this will be the marking score. (can then try to weight this later)

    therefore will return a list of 'correctness' where each value in the list corresponds to one injected transit. 
    - some may be seen, others not. 

    input
    ----
    
    grp: array of all the information for all clssifications for one tic_id


    '''
    # see if the TIC ID has a planet by seeing if it's in the ground truths
    # information needed from whole grp: min/max pix and time, TIC. (TIC_ID, min_pix, max_pix and min_time, max_time)
    # individual informatio needed: location and width of marked transit (xvals, width and frame)


    tic = int(grp['TIC_ID'].iloc[0])
    min_time = float(grp['min_time'].iloc[0])
    max_time = float(grp['max_time'].iloc[0])

    # information that is the same for all classifications in this group yay 

    # a list of the transits for that tic and a list of a single value of the width of all those transits.
    transits_time, widths  = transit_loc_per_tic(tic, min_time,max_time) # INJECTED (true transits)
    
    # if tic not inlist (no transit), returns [0,0], [0]
    if widths == [0]: # if rhere is no transit, return 999.9

        number_marked_false = 0 

        for row_n, row in grp.iterrows():
            try:
                number_marked_false += len(row['xvals'])
            except:
                number_marked_false += 0 

        return False, number_marked_false
        
        # time is a list of the times of the injected transits
    else:    
        # information that is the same for all classifications in this group yay 

        min_pix = float(grp['min_pix'].iloc[0])
        max_pix = float(grp['max_pix'].iloc[0])
        
        transits_pixels = time_to_pixel(transits_time, min_time, max_time) # get a list back

        transits_widths = deltat_to_pixel(widths, min_time, max_time) # get a list back

        width_pixel = transits_widths

        # make a list of zeros of length of the number of injected transit - 
        # true or false because this list will be updated on whether people fidm transits or not

        frac_marked_true = np.zeros_like(transits_pixels)
        number_marked_false = 0  # the number of transits marked where there were no transits.

        for row_n, row in grp.iterrows(): #iterate through all the rows in the group
            
            xvals = row['xvals']
            widths = row['width']

            if xvals != None:

                # need to loop through the marked transits and the injected transits (both boxes)
                # compare the edges of the boxes 
                # if marked min or max is within injected min or max
                # we are counting, the correctness of each injected transit        
                # loop through the known injected transits

                for i,inj in enumerate(transits_pixels):
                    
                    overlap = 0

                    for m, marked in enumerate(xvals):

                        # if they just drew a massive box, don't count it - 
                        # everything is in pixels therfore make as a functiom of the size of the image

                        if widths[m] < (max_pix - min_pix)/1:
                        
                             ## how many transits were marked in places where there were no transits?
                             # if the marked transit overlaps with any of the  injected ones, add one to  overlap.
                             # if the marked transit

                            # how many marked transits overlap with the injected transits?
                            # if satisfied break out of the loop 
                            if (inj + (width_pixel[0])/2) > (marked - widths[m]/2) > (inj - (width_pixel[0])/2) or (inj + (width_pixel[0])/2) > (marked + widths[m]/2) > (inj - (width_pixel[0])/2):
                                frac_marked_true[i] += 1
                                overlap += 1
                                break # have thae break condition so that each marking can only be counted as correct once
                

                    if overlap == 0:  # if there were no overlaps, then the transit was marked falsely
                        number_marked_false += 1

        frac_marked_true = [i/len(grp) for i in frac_marked_true] # calculate the numbers marked as a fraction of the total people who saw that lightcurve
        # this is fraction of everyonem should be only out of people who said yes? 


        return frac_marked_true, number_marked_false



