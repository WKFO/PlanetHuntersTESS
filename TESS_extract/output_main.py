#Python 2.7.9 (default, Apr  5 2015, 22:21:35)
# full env in environment.yml
from __future__ import print_function, division, absolute_import

import sys, os
import ujson
import random
import gc
from argparse import ArgumentParser

from astropy.table import Table
import numpy as np  # using 1.10.1
import pandas as pd  # using 0.13.1
#import datetime
#import dateutil.parser
import scipy
from PyAstronomy.pyTiming import pyPDM  # for the pdm analysis
import pwkit
import matplotlib.pyplot as plt

#import vote_distributions
import aggregate_class
import weighting
import get_info
import transit_gauss

# input into the command line
# python output_main.py


'''
This is a full aggregation of the Rel01 project, including user weighting.


This work flow has chopped each 27.8 day lightcurve into 4 lightcurves to make 4 inficidual subjects. 

It starts off with a simple Yes/No question. The weighting is based on this. 
It then goes on to look at the markign of the transits - each tarnsit is marked with a column tool. 
Each column is modelled as a gaussian curve with the width being the FWHM. 
'''

# file with raw classifications (csv) needed
# put this way up here so if there are no inputs we exit quickly before even trying to load everything else



try:
    classfile_in = "/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/8349/planet-hunters-tess-classifications.csv" # input file
except:
    #classfile_in = 'exoplanet-explorers-classifications.csv'

    print("      classifications_infile is a Zooniverse (Panoptes) classifications data export CSV.")
    print("      weight_class is 1 if you want to calculate and apply user weightings, 0 otherwise.")
    print("      aggregations_outfile is the name of the file you want written. If you don't specify,")
    print("      the filename is %s by default." % outfile_default)
    sys.exit(0)


###################### Define files and settings first ###################### 
#############################################################################

# define the active workflow - we will ignore all classifications not in this workflow

active_workflow_id = 8349
active_workflow_major = 12
apply_weight = 3

# what data release are we looking at e.g. ETE-6, Rel00, Rel01
data_release = "Rel01"
rankfile_stem = 'subjects_ranked_by_weighted_'


# do we want sum(weighted vote count) = sum(raw vote count)?
normalise_weights = True
# do we want to write an extra file with just classification counts and usernames
# (and a random color column, for treemaps)?
counts_out = True

# this is hard coded for the test verions - in later versions will only have one WF so will not need this. 

counts_out_file = '/Users/Nora/Documents/research/TESS/planethunters/{}/output/outfiles/class_counts_{}_{}.csv'.format(data_release,active_workflow_id, active_workflow_major)
outfile = '/Users/Nora/Documents/research/TESS/planethunters/{}/output/outfiles/planet_aggregations_{}_{}.csv'.format(data_release,active_workflow_id, active_workflow_major)



#################################################################################


# Print out the input parameters just as a sanity check
print("Computing aggregations using:")
print("infile: %s" % classfile_in)
print("weighted? %d" % apply_weight)
print("Will print to %s after processing." % outfile)


#********************************************************************************
#################################################################################
#********************************************************************************

                            # Begin the main work

#********************************************************************************
#################################################################################
#********************************************************************************

# LOAD THE DATA
# print("Reading classifications from %s ..." % classfile_in)
classifications = pd.read_csv(classfile_in) # this step can take a few minutes for a big file


# first, extract the started_at and finished_at from the annotations column
# print("Making new columns and getting user labels...")
classifications['metadata'] = [ujson.loads(q) for q in classifications.metadata]

classifications['subject_dimensions'] = [q['subject_dimensions'] for q in classifications.metadata]

# print(classifications['subject_dimensions'])
classifications['dimension_pixel'] = [get_info.get_dimension_pix(q) for q in classifications['subject_dimensions']]  # get the true/false for planet+EB and stars


#----------------------------------------------------------------------------
# Filter out classifications by me... (I don't do them properly when testing, sometimes neither does Grant so maybe filter him out too)
#----------------------------------------------------------------------------

#classifications = classifications[classifications.user_name != "nora.eisner"]
#classifications = classifications[classifications.user_name != "mrniaboc"]    # Grant 

# extract the classification annotations
classifications['annotation_json'] = [ujson.loads(q) for q in classifications.annotations]

# get the classifications of all of them - this is a list if they classified (yes), and a an empty list if they didn't (no)
classifications['planet_classification0'] = [q[0]['value'] for q in classifications.annotation_json]  

#make a new column where None classifications are calle NoClass
classifications['planet_classification'] = classifications.apply(get_info.yes_no_none, axis=1)

#delete the old planet_classification0 column 
del classifications['planet_classification0']


#-------------------------------------------------------------

#-------------------------------------------------------------
# Only use LIVE classifications  --- EVERYTHING IS NOT LIVE AT THE MOMENT
#-------------------------------------------------------------

# would that we could just do q['live_project'] but if that tag is missing for
# any classifications (which it is in some cases) it crashes

#classifications['live_project']  = [get_info.get_live_project(q) for q in classifications.metadata]
#
## if this line gives you an error you've read in this boolean as a string
## so need to convert "True" --> True and "False" --> False
#class_live = classifications[classifications.live_project].copy()
#n_class_thiswf = len(classifications)
#n_live = sum(classifications.live_project)
#n_notlive = n_class_thiswf - n_live
#print(" Removing %d non-live classifications..." % n_notlive)
#
## don't make a slice but also save memory
#classifications = pd.DataFrame(class_live)  # new panda dataframe with only the live classifications
#del class_live
#gc.collect()  # delete the deleted stuff (i.e. empty trash)

#-------------------------------------------------------------


# discard classifications not in the active workflow  
#-------------------------------------------------------------

print("Picking classifications from the active workflow (id %d, version %d.*)" % (active_workflow_id, active_workflow_major))
# use any workflow consistent with this major version, e.g. 6.12 and 6.23 are both 6 so they're both ok
# also check it's the correct workflow id
the_active_workflow = [int(q) == active_workflow_major for q in classifications.workflow_version]
this_workflow = classifications.workflow_id == active_workflow_id
in_workflow = this_workflow & the_active_workflow
# note I haven't saved the full DF anywhere because of memory reasons, so if you're debugging:
# classifications_all = classifications.copy()
classifications = classifications[in_workflow]  # the new classification now only has the dedicated WF and version

#-------------------------------------------------------------
print ("Number of yes' and no's...")
print (classifications['planet_classification'].value_counts())
#-------------------------------------------------------------


# get uniqueness of classifications
#----------------------------------

classifications['started_at_str']  = [q['started_at']  for q in classifications.metadata]
classifications['finished_at_str'] = [q['finished_at'] for q in classifications.metadata]

classifications['time_spent'] = [get_info.total_time(q) for q in classifications['started_at_str finished_at_str'.split()].iterrows()]

# print the mean, median, total time spent on a classification
print ("mean time spent on clssificaton: {}".format(np.mean(classifications['time_spent'][classifications['time_spent'] < 10])))
print ("median time spent on clssificaton: {}".format(np.median(classifications['time_spent'][classifications['time_spent'] < 10])))
print ("total time spent on clssificaton: {}".format(np.sum(classifications['time_spent'][classifications['time_spent'] < 10])))

# we need to set up a new user id column that's login name if the classification is while logged in,
# session if not (right now "user_name" is login name or hashed IP and, well, read on...)
classifications['user_label'] = [get_info.get_alternate_sessioninfo(q) for q in classifications['user_name metadata'.split()].iterrows()]
classifications['created_day'] = [q[:10] for q in classifications.created_at]

# Get subject info into a format we can actually use
classifications['subject_json'] = [ujson.loads(q) for q in classifications.subject_data]


#-----------------------------------
# Classification yes/no information   
#-----------------------------------

print("Getting classification info...")

# Get annotation info into a format we can actually use

# marked transits == there is a transit. 
# nothing marked == there are no transits.

#-----------------------------------
# create a weight parameter but set it to 1.0 for all classifications (unweighted) - may change later
classifications['weight'] = [1.0 for q in classifications.workflow_version]
# also create a count parameter, because at the time of writing this .aggregate('count') was sometimes off by 1
classifications['count'] = [1 for q in classifications.workflow_version]

# ...we will use these later
#-----------------------------------


print("Extracting filenames and subject types...")

'''
 extract whether this is a known planet or a candidate that needs classifying - that info is in the "planet" column in the subject metadata (where the ? can be either "#", "!" or "//").
 do this after you choose a workflow because #Class doesn't exist for the early subjects so it will break
 also don't send the entirety of classifications into the function, to save memory
'''

# IS THERE A PLANET OR AN EB IN THE SIMULATED DATA - alter deoending on whether we want EB to count as 'yes found'
classifications['subject_type'] = [get_info.get_subject_type(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]  # get the true/false for planet+EB and unknown
classifications['subject'] = [get_info.get_subject(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]            # get the actual subject

classifications['candidate'] = [get_info.get_filename(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]
classifications['TIC_ID'] = [get_info.get_subject_TIC(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]

# Extract the pixels and times of the frames:
classifications['min_pix'] = [get_info.get_min_pix(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]
classifications['max_pix'] = [get_info.get_max_pix(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]

classifications['min_time'] = [get_info.get_min_time(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]
classifications['max_time'] = [get_info.get_max_time(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]


# REMEMBER TO TURN THIS ON!!!!!!!
#-----------------------------------
#-----------------------------------
#-----------------------------------

# WHEN LIVE, THIS WILL WORK - check what the correct 'true' will be
# True = is it a simulated LC, False = not simulated
classifications['sim'] = [get_info.get_sim(q) for q in classifications['subject_ids subject_json'.split()].iterrows()]

#-----------------------------------
#-----------------------------------
#-----------------------------------


# import the ground truth data to get the SNR for the planets. 
#-----------------------------------


#----------------
# Assess the SNR   # for the simulated LCs if not SIM data, then return None
#----------------

print("Getting ground truth information for SNR values...")
# this function gets the SNR from the ground data for each planet (there is no SNR value for not-planets).

classifications_snr=classifications[classifications.subject_type ==True] # True mans that it's a simulated data

SNR_DF = get_info.get_SNR(classifications_snr)

# add the SNR datafile to the main one by merging the two data frames
classifications_old = classifications.copy()
classifications = classifications_old.merge(SNR_DF,how='outer', left_index=True, right_index=True)



#---------------------------------
# Location of the Marked Transits
#---------------------------------

print("Finding the location of the marked transits...")

# extract location, width and frame (for zoomed version) from the json
# these are all in the form of lists, as each user can mark more than one transit. 
# extract this information for each user and each classification (no groupping yet)
#x_DF, width_DF, frame_DF = get_info.get_marking(classifications)

marking_DF = get_info.get_marking(classifications)


# add two new columns to the DF: width and cenre
classifications_old = classifications.copy()
classifications = classifications_old.merge(marking_DF,how='outer', left_index=True, right_index=True)


#---------------------------------------------------

# this might be useful for a sanity check later
# first_class_day = np.min(classifications.created_day).replace(' ', '')
# last_class_day  = np.max(classifications.created_day).replace(' ', '')
# for some reason this is reporting last-classification dates that are days after the actual last
# classification. Might be because this is a front-end reporting, so if someone has set
# their computer's time wrong we could get the wrong time here.
# could fix that by using created_at...
# last_class_time = np.max(classifications.finished_at_str)[:16].replace(' ', '_').replace('T', '_').replace(':', 'h')+"m"

last_class_time = np.max(classifications.created_at)[:16].replace(' ', '_').replace('T', '_').replace(':', 'h')+"m"



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
#######################################################
#         Apply weighting function (or don't)         #
#######################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #


classifications['seed'] = [0 for q in classifications.weight]
classifications['is_gs'] = [0 for q in classifications.weight]

if apply_weight > 0:

    print("  Computing user weights...")

    # if it breaks on this and gives a "str + int" error check these are
    # actually booleans and not just strings that say "True" or "False"
    # and note setting an array of those strings with a .astype(bool)

    # this is a list of TRUE = transit, FALSE = no transit





    # CHANGE LATER !!!!!!!!!!!!!!!!!!!!!!! TO SIM
    is_known = classifications.subject_type.values  
    #is_known = classifications.sim.values  





    #is_candidate = np.invert(is_known)
    # if it's a non-gold-standard classification, mark it
    classifications.loc[is_known, 'is_gs'] = 1


    if apply_weight < 2:
        # the simple case of incrementing or decrementing the seed weight
        true_incr   = 1.0  # upweight if correct
        # downweight only a tiny (constant) amount if users missed a sim
        false_incr = -0.01
    else:
        # this will be scaled by how easy the sim is to spot
        # so this is the max upweight if ~all users missed it and one spotted it
        true_incr   = 2.0
        # this is the max downweight if ~all users spotted it and one missed it
        false_incr = -1.0

        '''
          I am really not convinced this is the right weighting scheme; it's probably worth testing various things out, including trying to make an actual confusion matrix.
          Well, it would be, if we had true negative sims. But we don't.
        '''

    # find the correct classifications of known planets
    TP  = (is_known) & (classifications.subject_type == True) & (classifications.planet_classification == 'Yes')
    
    # find the incorrect classifications of known planets   --- this currently includes None 
    FN = (is_known) & (classifications.subject_type == True) & (classifications.planet_classification == 'No')

    # if we're weighting differently based on how hard each particular sim is to spot
    # (measured by what fraction of people spotted it)
    # then we need to aggregate just the sims to get the unweighted fractions
    # and then update the seeds for those
    
    if apply_weight == 3:
        print("Weighting of sims based on sim difficulty (this takes longer)...")
        by_subj_sims = classifications[is_known].groupby('candidate') # if it's known to be a planet, group by candidate

        sims_agg = by_subj_sims['weight count planet_classification subject_type subject candidate TIC_ID SNR min_pix max_pix min_time max_time'.split()].apply(aggregate_class.aggregate_class)

        # replace NaN with 0.0
        sims_agg.fillna(value=0.0, inplace=True)

        '''
        We want to adjust the seed for each subject's weight based on what
        fraction of classifiers spotted the planet in the sim subject
        if we're weighting Alice's classification and Alice said Yes to a sim
        (i.e. got it right), but everyone else got it right too, then this
        subject doesn't give us much information about her ability, so the
        seed should be low. On the other hand, if she spotted it when most
        others missed it, she should get a nice upweight as a reward,
        so the multiplier for the ok_class group is actually p_No.
        likewise, if Bob failed to spot a sim but so did everyone else, then
        the sim is just impossible to spot and Bob shouldn't be penalized for
        that, whereas if everyone spotted it but Bob, Bob should be downweighted
        more heavily; so the multiplier for the oops_class group is p_Yes.
        We don't have true-negative sims, so if someone correctly spots a very
        difficult sim we can't tell whether they are really that astute or
        whether they're just really optimistic about transits.
        So this weighting scheme has the effect of rewarding the optimists with higher weightings, and I'm unconvinced that's the right way to go. After
        some spot checking I think it's at least possible we should be trying to
        get the false positive rate down a bit.
        '''

        sims_agg['fseed_right']   = sims_agg['p_No']
        sims_agg['fseed_false']   = sims_agg['p_Yes']

        # now we have multipliers for each sim subject, add them back to the classifications dataframe

        # add two new columns into the data frame that look at the yes and No probabilities f
        classifications_old = classifications.copy()
        classifications = pd.merge(classifications_old, sims_agg['fseed_right fseed_false'.split()], how='left', left_on='candidate', right_index=True, sort=False, suffixes=('_2', ''), copy=True)

        del classifications_old

        #classifications['fseed_right fseed_false'.split()].fillna(value=0.0, inplace=True)
        # set the individual seeds
        
        # changing the seeds for each ok_class. Classification is correct then multiple the percntage of false positive with the ok incr.
        
        # where is known and correct , 'seed' == the incrememend (positive) times the percentage of people who got it right
        classifications.loc[TP, 'seed'] = true_incr * classifications['fseed_right'][TP] 

        # where is known and incorect , 'seed' == the incremement (negative) times the percentage of people who got it wrong 
        classifications.loc[FN, 'seed'] = false_incr * classifications['fseed_false'][FN]
        
        # by this point we have generated a seed FOR EACH CLASSIFICATION, not each user.
        # this seed is based on whether their response agrees with the overall concensus
        # more positive if it's correct but lots of people got it wrong
        # more negative if they got it wrong but everyone else got it right

    # then group classifications by user name, which will weight logged in as well as not-logged-in (the latter by session)
    by_user = classifications.groupby('user_label')

    # For each user, sum all of their seeds for each classification that they did. 
    # This goes into the exponent for the weight
    # This is their overall seed

    user_exp = by_user.seed.aggregate('sum')

    # then set up the DF that will contain the weights etc, and fill it
    user_weights = pd.DataFrame(user_exp)
    user_weights.columns = ['seed'] # label the summed seeds "seed"

    user_weights['user_label'] = user_weights.index
    user_weights['user_id'] = by_user['user_id'].head(1)
    user_weights['nclass_user'] = by_user['count'].aggregate('sum')
    user_weights['n_gs'] = by_user['is_gs'].aggregate('sum')
    user_weights['weight'] = [weighting.assign_weight(q, apply_weight) for q in user_weights.iterrows()]

    # user_weights['weight'] = [assign_weight_old(q) for q in user_exp]

    # if you want sum(unweighted classification count) == sum(weighted classification count), do this
    if normalise_weights:
        user_weights['weight_unnorm'] = user_weights['weight'].copy()
        user_weights.weight *= float(len(classifications))/float(sum(user_weights.weight * user_weights.nclass_user))

    # weights are assigned, now need to match them up to the main classifications table
    # making sure that this weight keeps the name 'weight' and the other gets renamed (suffixes flag)
    # if assign_weight == 0 then we won't enter this loop and the old "weights" will stay
    # as they are, i.e. == 1 uniformly.

    classifications_old = classifications.copy()
    classifications = pd.merge(classifications_old, user_weights, how='left',on='user_label', sort=False, suffixes=('_2', ''), copy=True)
    del classifications_old
        

else:
    # just make a collated classification count array so we can print it to the screen
    by_user = classifications.groupby('user_label')
    user_exp = by_user.seed.aggregate('sum')
    user_weights = pd.DataFrame(user_exp)
    user_weights.columns = ['seed']
    #user_weights['user_label'] = user_weights.index
    user_weights['nclass_user'] = by_user['count'].aggregate('sum')
    user_weights['n_gs'] = by_user['is_gs'].aggregate('sum')
    
    # UNWEIGHTED
    user_weights['weight'] = [1 for q in user_exp]


gc.collect()


# grab basic stats
n_subj_tot  = len(classifications.subject_data.unique())
by_subject = classifications.groupby('candidate')
subj_class = by_subject.created_at.aggregate('count')
all_users  = classifications.user_label.unique()
n_user_tot = len(all_users)
n_user_unreg = sum([q.startswith('not-logged-in-') for q in all_users])
last_class_id = max(classifications.classification_id)

# obviously if we didn't weight then we don't need to get stats on weights
if apply_weight > 0:
    user_weight_mean   = np.mean(user_weights.weight)
    user_weight_median = np.median(user_weights.weight)
    user_weight_25pct  = np.percentile(user_weights.weight, 25)
    user_weight_75pct  = np.percentile(user_weights.weight, 75)
    user_weight_min    = min(user_weights.weight)
    user_weight_max    = max(user_weights.weight)

nclass_mean   = np.mean(user_weights.nclass_user)
nclass_median = np.median(user_weights.nclass_user)
nclass_tot    = len(classifications)

user_weights.sort_values(['nclass_user'], ascending=False, inplace=True)


#------------------------------------------------------
#                  Transit Marking                    #
#------------------------------------------------------

# this cannot be done with TIC_IDs for this project - needs to be done with the image name!

# evaluate transit markings as a gaussian
#-----------------------------------------

candidates = classifications['candidate'].unique()
data = []

for candidate in candidates:
    classifications_of_candidate = classifications[classifications['candidate'] == candidate]
    x_eval, candidate_pdf = transit_gauss.xval_hist_gauss(classifications_of_candidate)
    x_peaks, pdf_peaks = transit_gauss.get_distribution_peaks(x_eval, candidate_pdf)
    
    x_peak_max, pdf_peak_max = transit_gauss.get_max_peak(x_eval, candidate_pdf) # the highest peak - can use this for a ranking as an alternative to yes/no

    #get_info.DB_scan(classifications_of_candidate)

    # return either 999.9 if mo transits exist, or an array of the number of correct markings per transit.

    frac_marked_true,number_marked_false = transit_gauss.peaks_right_wrong(classifications_of_candidate)

    data.append({'candidate': candidate, 'x_eval': list(x_eval), 'pdf': list(candidate_pdf), 'x_peaks':list(x_peaks), 'pdf_peaks':list(pdf_peaks), 'x_peak_max':x_peak_max, 'pdf_peak_max':pdf_peak_max,'frac_marked_correct':frac_marked_true, 'number_marked_false':number_marked_false})

peak_df = pd.DataFrame(data)


# evaluate transit markings as a boxes
#-----------------------------------------

candidates = classifications['candidate'].unique()
data_box = []

for candidate in candidates:
    classifications_of_candidate = classifications[classifications['candidate'] == candidate]
    x_eval_box, candidate_pdf_box = transit_gauss.xval_hist_box(classifications_of_candidate)
    x_peaks_box, pdf_peaks_box = transit_gauss.get_distribution_peaks(x_eval_box, candidate_pdf_box)
    
    x_peak_max_box, pdf_peak_max_box = transit_gauss.get_max_peak(x_eval_box, candidate_pdf_box) # the highest peak - can use this for a ranking as an alternative to yes/no

    #get_info.DB_scan(classifications_of_candidate)

    # return either 999.9 if mo transits exist, or an array of the number of correct markings per transit.
    data_box.append({'candidate': candidate, 'x_eval_box': list(x_eval_box), 'pdf_box': list(candidate_pdf_box), 'x_peaks_box':list(x_peaks_box), 'pdf_peaks_box':list(pdf_peaks_box), 'x_peak_max_box':x_peak_max_box, 'pdf_peak_max_box':pdf_peak_max_box})

peak_df_box = pd.DataFrame(data_box)


# use A DB scan clustering algorithm to get centres and counts, and weighted counts
#-----------------------------------------

data_db = []
for candidate in candidates:
    classifications_of_candidate = classifications[classifications['candidate'] == candidate]
    db_peak,db_count,db_count_weighted = get_info.DB_scan(classifications_of_candidate)
    data_db.append({'candidate': candidate,'db_peak':db_peak,'db_count':db_count, 'db_count_weighted':db_count_weighted})

db_df = pd.DataFrame(data_db)

# try only a couple of DB scans:
#for candidate in candidates[20:22]:
#
#    classifications_of_candidate = classifications[classifications['candidate'] == candidate]
#    get_info.DB_scan(classifications_of_candidate)



#------------------------------------------------------

## ## ## ## ## ## ##             ## ## ## ## ## ## ## #
#######################################################
#            Print out basic project info             #
#######################################################
## ## ## ## ## ## ##             ## ## ## ## ## ## ## #


print("%d classifications from %d users, %d registered and %d unregistered.\n" % (nclass_tot, n_user_tot, n_user_tot - n_user_unreg, n_user_unreg))
print("Mean n_class per user %.1f, median %.1f." % (nclass_mean, nclass_median))
if apply_weight > 0:
    print("Mean user weight %.3f, median %.3f, with the middle 50 percent of users between %.3f and %.3f." % (user_weight_mean, user_weight_median, user_weight_25pct, user_weight_75pct))
    print("The min user weight is %.3f and the max user weight is %.3f.\n" % (user_weight_min, user_weight_max))
    cols_print = 'nclass_user weight'.split()
else:
    cols_print = 'nclass_user'
# don't make this leaderboard public unless you want to gamify your users in ways we already know
# have unintended and sometimes negative consequences. This is just for your information.
print("Gini coefficient for project: %.3f" % weighting.gini(user_weights['nclass_user']))

# If you want to print out a file of classification counts per user, with colors
# for making a treemap
# honestly I'm not sure why you wouldn't want to print this, as it's very little
# extra effort AND it's the only place we record the user weights



#######################################################
#                   Write to files                    #
#######################################################


if counts_out == True:
    print("Printing classification counts to %s..." % counts_out_file)
    user_weights['color'] = [get_info.randcolor(q) for q in user_weights.index]
    user_weights.to_csv(counts_out_file)


#------------------------------------------------------
# Aggregate classifications, unweighted and weighted  #
#------------------------------------------------------

print("\nAggregating classifications...\n")

class_agg = by_subject['weight count planet_classification subject_type subject candidate TIC_ID SNR min_pix max_pix min_time max_time'.split()].apply(aggregate_class.aggregate_class)
# if p_Ans = 0.0 the function skipped those values, so replace all the NaNs with 0.0
class_agg.fillna(value=0.0, inplace=True)


# class_agg is grouped by TIC
# add the peak/position information to the main data frame which is stored to the file for further evaluation.
# remember that these peak positions are still in terms of pixels - will change to be in terms of time in the next step. 
class_agg_old = class_agg.copy()
class_agg = class_agg_old.merge(peak_df,on='candidate')

class_agg_old = class_agg.copy()
class_agg = class_agg_old.merge(peak_df_box,on='candidate')

class_agg_old = class_agg.copy()
class_agg = class_agg_old.merge(db_df,on='candidate')


# add value-added columns
# let people look up the subject on Talk directly from the aggregated file

class_agg['link'] = ['https://www.zooniverse.org/projects/ianc2/exoplanet-explorers/talk/subjects/'+str(q) for q in class_agg.index]
# and provide the exofop link too
# the int(float(q)) thing is in case it reads as a string, which breaks int(q)
# class_agg['exofop'] = ['https://exofop.ipac.caltech.edu/k2/edit_target.php?id=%d' % int(float(q)) for q in class_agg.candidate]
# after we do the merges below the new indices might not be linked to the subject id, so save it explicitly
class_agg['subject_ids'] = [str(q) for q in class_agg.index]

# make the list ranked by p_Yes_weight
class_agg.sort_values(['subject_type','p_Yes_weight'], ascending=False, inplace=True)

print("Writing aggregated output to file %s...\n" % outfile)
pd.DataFrame(class_agg).to_csv(outfile)


#######################################################

# Now make files ranked by p_Yes, one with all subjects classified and one with only candidates

# /Users/vrooje/anaconda/bin/ipython:1: FutureWarning: sort(columns=....) is deprecated, use sort_values(by=.....)
# !/bin/bash /Users/vrooje/anaconda/bin/python.app
#class_agg.sort('p_Yes_weight', ascending=False, inplace=True)
class_agg.sort_values(['p_Yes_weight'], ascending=False, inplace=True)

# I'd rather note the last classification date than the date we happen to produce the file
# rightnow = datetime.datetime.now().strftime('%Y-%M-%D_%H:%M')
# rankfile_all = rankfile_stem + rightnow + ".csv"
rankfile_all = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/all_%s_%s_lastid_%d.csv' % (rankfile_stem, last_class_time, last_class_id)


# there go those hard-coded columns again
#rank_cols = ['subject_ids', 'candidate', 'p_Yes_weight', 'count_weighted', 'p_Yes', 'count_unweighted', 'subject_type', 'link', 'user_tag', 'HTRU-N File']
rank_cols = ['subject_ids', 'TIC_ID', 'candidate', 'p_Yes_weight', 'count_weighted', 'p_Yes', 'count_unweighted', 'subject_type', 'subject', 'link', 'SNR', ]


print("Writing full ranked list to file %s...\n" % rankfile_all)
# write just the weighted yes percentage, the weighted count, the subject type, and the link to the subject page
# the subject ID is the index so it will be written anyway
#pd.DataFrame(class_agg[rank_cols]).to_csv(rankfile_all)


rankfile = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/outfiles/nonsim_allsubj_%s_%s_lastid_%d.csv' % (rankfile_stem, last_class_time, last_class_id)
print("Writing candidate-only ranked list to file %s...\n" % rankfile)
# also only include entries where there were at least 5 weighted votes tallied
# and only "cand" subject_type objects
classified_candidate = (class_agg.count_weighted > 5) & (class_agg.subject_type == False)
#pd.DataFrame(class_agg[rank_cols][classified_candidate]).to_csv(rankfile)

#done.