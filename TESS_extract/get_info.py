from astropy.table import Table
import random
import pandas as pd  # using 0.13.1
from astropy.time import Time
# These are functions that extract information from the various JSONs that are
# included in the classification exports.
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt


def get_subject_type(q):

    try:
        if q[1].subject_json[str(q[1].subject_ids)]['#Target'] == 'planet' or q[1].subject_json[str(q[1].subject_ids)]['#Target'] == 'EB':
            return True
        else:
            return False
    except:
        if q[1].subject_json[str(q[1].subject_ids)]['!Target'] == 'planet' or q[1].subject_json[str(q[1].subject_ids)]['!Target'] == 'EB':
            return True
        else:
            return False


def get_sim(q):

        #if q[1].subject_json[str(q[1].subject_ids)]['#sim'] == 'True':
        #    return True
        #else:
        #    return False
        try:
            if q[1].subject_json[str(q[1].subject_ids)]['#sim'] == 'True':
                return True
            else:
                return False
        except:
            return False


def get_subject(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['#Target']
    except:
        return q[1].subject_json[str(q[1].subject_ids)]['!Target']

def get_subject_TIC(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['!TIC ID'] 
    except:
        return q[1].subject_json[str(q[1].subject_ids)]['TIC ID'] 


def get_min_pix(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['#xmin_pix'] 
    except:
        return q[1].subject_json[str(q[1].subject_ids)]['#xmin_pix_main'] 

def get_max_pix(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['#xmax_pix']
    except:
        return q[1].subject_json[str(q[1].subject_ids)]['#xmax_pix_main']

def get_min_time(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['#xmin_time']
    except: 
        return q[1].subject_json[str(q[1].subject_ids)]['#xmin_time_main'] 

def get_max_time(q):
    try:
        return q[1].subject_json[str(q[1].subject_ids)]['xmax_time'] 
    except:
        try:
            return q[1].subject_json[str(q[1].subject_ids)]['#xmax_time_main'] 
        except:
            try:
                return q[1].subject_json[str(q[1].subject_ids)]['#xmax_time'] 
            except:
                print("move on")


def get_dimension_pix(q):

    try:
        return q['clientWidth'] 
    except:
        return 999.999


def get_filename(q):
    try:
        try:
            return q[1].subject_json[str(q[1].subject_ids)]['#Image']
        except:
            return q[1].subject_json[str(q[1].subject_ids)]['#Image0']
    except:
        q[1].subject_json[str(q[1].subject_ids)]['!Image0']


# get number of gold-standard classifications completed by a user (used if weighting)
def get_n_gs(thegrp):
    return sum(pd.DataFrame(thegrp).seed != 0)


# Something went weird with IP addresses, so use more info to determine unique users
# Note the user_name still has the IP address in it if the user is not logged in;
# it's just that for this specific project it's not that informative.
# Note: the above was for Pulsar Hunters and IPs worked fine in Exoplanet Explorers
# so we don't really need this, but keep it in case you ever want to try to identify
# unique users based on something other than IP address


def get_alternate_sessioninfo(row):

    # if they're logged in, save yourself all this trouble
    if not row[1]['user_name'].startswith('not-logged-in'):
        return row[1]['user_name']
    else:
        metadata = row[1]['metadata']
        # IP + session, if it exists
        # (IP, agent, viewport_width, viewport_height) if session doesn't exist
        try:
            # start with "not-logged-in" so stuff later doesn't break
            return str(row[1]['user_name']) +"_"+ str(metadata['session'])
        except:
            try:
                viewport = str(metadata['viewport'])
            except:
                viewport = "NoViewport"

            try:
                user_agent = str(metadata['user_agent'])
            except:
                user_agent = "NoUserAgent"

            try:
                user_ip = str(row[1]['user_name'])
            except:
                user_ip = "NoUserIP"

            thesession = user_ip + user_agent + viewport
            return thesession


def total_time(row):

    # required format: 2017-10-27 23:03:09
    started_t = (row[1]['started_at_str'][0:10]) + ' ' + (row[1]['started_at_str'][11:-2])
    finished_t = (row[1]['finished_at_str'][0:10]) + ' ' + (row[1]['finished_at_str'][11:-2])

    t_start = Time(started_t, format='iso', scale='utc') # convert to JD
    t_finish = Time(finished_t, format='iso', scale='utc') # convert to JD

    
    start_JD = t_start.jd
    finish_JD = t_finish.jd

    delta_t = (finish_JD - start_JD) * 24 * 60 # get into minutes

    return delta_t



def get_live_project(metadata):
    try:
        return metadata['live_project']
    except:
        # apparently some subject metadata doesn't have this? dunno?
        return False


import numpy as np
from scipy.stats import norm


def get_net_voting_disget_max_peaktribution(xvals, widths, weight, x_eval):
    """
    Find smooth pdf of contribution of user

    Args:
    ----
        xvals (list): x's of transits identified by user
        widths (list): widths of transits identified by user
        weight (float): relative weight of user, for Gaussian area

    Returns:
    -------
        (np.array): net voting distribution, evaluated at x_eval points
    """

    assert len(xvals) == len(widths)

    try:

        classification_pdf = np.zeros_like(x_eval)

        for classification_n in range(len(xvals)):
            x = xvals[classification_n]
            width = widths[classification_n]
            classification_pdf += norm.pdf(x_eval, loc=x, scale=width/2.355)   # divide the width by two to use the FWHM instead of the standdard deviation
        
        return classification_pdf * weight  # scale up by user  #user_vote_pdf
    
    except:
        return np.zeros_like(x_eval)
        

def get_net_voting_distribution_box(xvals, widths, weight, x_eval):
    """
    Find boxes

    Args:
    ----
        xvals (list): x's of transits identified by user
        widths (list): widths of transits identified by user
        weight (float): relative weight of user, for Gaussian area

    Returns:
    -------
        (np.array): net voting distribution, evaluated at x_eval points
    """


    assert len(xvals) == len(widths)


    try:

        classification_pdf = np.zeros_like(x_eval)

        for classification_n in range(len(xvals)):  # list of markings from one person 
            x = xvals[classification_n]  # centre of box
            width = widths[classification_n]  # width of box
            height = weight

            for i in range(int(x-(width/2)),int(x+(width/2))):

                classification_pdf[i] =+ height

        return classification_pdf  # scale up by user  #user_vote_pdf
    
    except:
        return np.zeros_like(x_eval)
        


def get_SNR(classifications):

    pl = Table.read('/Users/Nora/Documents/research/TESS/ETE-6/ground_truth/ete6_planet_data_snr.txt', format='ascii',comment='#')

    pl_snr = pl['col16']
    pl_tic = pl['col1']

    snr = []
    idx = []

    pl_dict = {k: v for k, v in zip(pl_tic, pl_snr)}
    
    for i in classifications['TIC_ID']:

        if int(i) in list(pl_tic):
            snr.append(pl_dict[int(i)])
        else:
            snr.append(None)
    
    return pd.DataFrame(dict(SNR=snr),index=classifications['TIC_ID'].index.values)


def get_marking(classifications):

    """
    extract location, width and frame (for zoomed version) from the json
    these are all in the form of lists, as each user can mark more than one transit. 
    extract this information for each user and each classification (no groupping yet)

    arguments:
    ------
        classifcation (panda DF): DF of all the data (not grouped or aggregated)
    
    output:
    ------
        three panda data frames - for the width, xvals and which frame.

    """
    xvals = []
    width = []

    for k,cl in enumerate(classifications['annotation_json']): 
        
        if len(cl[0]['value']) > 0: # if someone marked something...
            
            x = []
            w = []
            f = []
            
            for i,v in enumerate(cl[0]['value']): 

                try:
                    x.append(v['x'])
                    w.append(v['width'])

                except:
                    x.append(None)
                    w.append(None)


            #x = [(i-45.89) * 0.039 for i in x]
            xvals.append(x)
            width.append(w)
               
        else:
            xvals.append(None)
            width.append(None)


    return pd.DataFrame(dict(xvals=xvals, width=width), index=classifications['TIC_ID'].index.values)

    #return pd.DataFrame(dict(xvals=xvals),index=classifications['TIC_ID'].index.values), pd.DataFrame(dict(width=width),index=classifications['TIC_ID'].index.values), pd.DataFrame(dict(frame=frame),index=classifications['TIC_ID'].index.values)


#################################################################################

# assign a color randomly if logged in, gray otherwise
# this is for making a treemap of your project's contributors
def randcolor(user_label):
    if user_label.startswith('not-logged-in-'):
        # keep it confined to grays, i.e. R=G=B and not too bright, not too dark
        g = random.randint(25,150)
        return '#%02X%02X%02X' % (g,g,g)
        #return '#555555'
    else:
        # the lambda makes this generate a new int every time it's called, so that
        # in general R != G != B below.
        r = lambda: random.randint(0,255)
        return '#%02X%02X%02X' % (r(),r(),r())


#################################################################################



def yes_no_none(row):

    """
    return whether people marked someting or not. 

    Pandas can't deal with our None values i.e. if people continue without classifying anything.
    We therefore need this function in order to make a new column that turns all None values into "NoClass" strings.
    """

    if len(row['planet_classification0']) == 0:
        val = "No"
    else:
        val = "Yes"
    
    return val


 
def DB_scan(grp):

    grp = grp.dropna(subset=['xvals', 'width']) # get rid of nan, None etc values (these will disrupt the code)

    min_pixels = 0
    max_pixels = 1193

    x_eval = np.linspace(min_pixels, max_pixels, 1193)
    grp_pdf = np.zeros_like(x_eval)
    
    xvals = [] # only care about the centre of the box if we are clustering
    weight = []

    for row_n, row in grp.iterrows():

        xvals_list = (row['xvals'])

        for x in xvals_list:
            xvals.append(x)
            weight.append((row['weight']))

    xvals_dic = [str(i) for i in xvals]
    dictionary = dict(zip(xvals_dic, weight))

    if len(xvals) != 0:
        yvals = [1] * len(xvals)
    
        X = np.column_stack((xvals, yvals))
    
        #X = StandardScaler().fit_transform(Xdb)
    
        db = DBSCAN(eps=20, min_samples=3).fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        
        # ------------------------------

        peak_position = []
        count = []
        count_weighted = []

        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]
        
            class_member_mask = (labels == k)
        
            xy = X[class_member_mask & core_samples_mask]
 
            if len(xy[:, 0]) > 0:
                peak_position.append(np.mean(xy[:, 0]))
                count.append(len(xy[:, 0]))

                count_w = 0
                for i in xy[:, 0]:
                    try:
                        w = dictionary['{}'.format(str(i))]
                        count_w += w
                    except:
                        w = dictionary['{}'.format(int(i))]
                        count_w += w              

                count_weighted.append(count_w)

        # get the count as a fraction of the total number of people who saw this LC
        count = [i/ len(grp['xvals']) for i in count]
        count_weighted = [i/ grp['weight'].sum() for i in count_weighted]

        if len(peak_position) == 0:
            return [0],[0],[0]
        else:
            return peak_position,count,count_weighted
    else:
        return [0],[0],[0]

def DB_scan_plot(grp):

    grp = grp.dropna(subset=['xvals', 'width']) # get rid of nan, None etc values (these will disrupt the code)

    min_pixels = 0
    max_pixels = 1193

    x_eval = np.linspace(min_pixels, max_pixels, 1193)
    grp_pdf = np.zeros_like(x_eval)
    
    xvals = [] # only care about the centre of the box if we are clustering
    weight = []

    for row_n, row in grp.iterrows():

        xvals_list = (row['xvals'])

        for x in xvals_list:
            xvals.append(x)
            weight.append((row['weight']))

    xvals_dic = [str(i) for i in xvals]
    dictionary = dict(zip(xvals_dic, weight))

    print (dictionary)

    if len(xvals) != 0:
        yvals = [1] * len(xvals)
    
        X = np.column_stack((xvals, yvals))
    
        #X = StandardScaler().fit_transform(Xdb)
    
        db = DBSCAN(eps=25, min_samples=3).fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        
        # ------------------------------

        peak_position = []
        count = []
        count_weighted = []

        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]
        
            class_member_mask = (labels == k)
        
            xy = X[class_member_mask & core_samples_mask]
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                     markeredgecolor='k', markersize=14)
            
            if len(xy[:, 0]) > 0:
                peak_position.append(np.mean(xy[:, 0]))
                count.append(len(xy[:, 0]))

                count_w = 0
                for i in xy[:, 0]:
                    print (i)
                    try:
                        w = dictionary['{}'.format(str(i))]
                        count_w += w
                    except:
                        w = dictionary['{}'.format(int(i))]
                        count_w += w              

                count_weighted.append(count_w)

            # outliers
            xy = X[class_member_mask & ~core_samples_mask]
            plt.xlim(0,1193)
            plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                     markeredgecolor='k', markersize=6)
        
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        plt.show()

