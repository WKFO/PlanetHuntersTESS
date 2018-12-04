import numpy as np


# The new weighting assignment function allows the user to choose between different weighting schemes
# though note the one in this function is not preferred for reasons explained below.
def assign_weight_old(seed):

    """
    Asess which weighting model is used. This changes the initial seed of each user. 

    Args:
    -----
        seed

    Returns:
    -----
        the intial seed

    """
    # keep the two seed cases separate because we might want to use a different base for each
    if seed < 0.:
        return np.max([0.05, pow(1.0025, seed)])
    elif seed > 0:
        return np.min([3.0, pow(1.0025, seed)])
    else:
        return 1.0


# assigns a weight based on a seed parameter
# The weight is assigned using the seed as an exponent and the number below as the base.
# The number is just slightly offset from 1 so that it takes many classifications for
# a user's potential weight to cap out at the max weight (3) or bottom out at the min (0.05).

# Note I'd rather this did a proper analysis with a confusion matrix etc but under a time crunch
# we went with something simpler.

def assign_weight(q, which_weight):

    """
    Asign a weight each user based on what the weighting system is.

    Args:
    -----
        q (panda DF): grouped by user - get one line per user. Analsyses one line at a time
        which_weight: which weighting system to use (3 takes into consideration how 'difficult' the transit is to see.)

    Returns:
    -----
        
    """

    # the floor weight for the case of which_weight == 2 or 3
    # i.e. someone who has seed = 0 will have this
    # seed = 0 could either be equal numbers right & wrong, OR that we don't have any information
    c0 = 0.5

    #user_weights['user_label'] = user_weights.index
    #user_weights['user_id'] = by_user['user_id'].head(1)
    #user_weights['nclass_user'] = by_user['count'].aggregate('sum')
    #user_weights['n_gs'] = by_user['is_gs'].aggregate('sum')

    seed         = q[1].seed   # summed seed for that person
    n_gs         = q[1].n_gs   # how many simulated lightcurves they have classified

    # Three possible weighting schemes:
    # which_weight == 1: w = 1.0025^(seed), bounded between 0.05 and 3.0
    # which_weight == 2: w = (1 + log n_gs)^(seed/n_gs), bounded between 0.05 and 3.0
    # which_weight == 3: same as == 2 but the seed is different; see below
    
    
    '''
     Weighting Scheme 1:
     this is an okay weighting scheme, but it doesn't account for the fact that someone might be prolific but not a very good classifier, and those classifiers shouldn't have a high weight.
     Example: Bob does 10000 gold-standard classifications and gets 5100 right, 4900 wrong. In this weighting scheme, Bob's weighting seed is +100, which means a weight of 1.0025^100 = 1.3, despite the fact that Bob's classifications are consistent with random within 1%.
     The weightings below this one would take the weight based on 100/10000, which is much better.
    '''
    if which_weight == 1:
        # keep the two seed cases separate because we might want to use a different base for each
        if seed < 0.:
            return np.max([0.05, pow(1.0025, seed)])
        elif seed > 0:
            return np.min([3.0, pow(1.0025, seed)])
        else:
            return 1.0

        '''
         Weighting Schemes 2 and 3:
         The passed seed is divided by the number of classifications.
         In weighting scheme 2, the seed has one fixed increment for correct classifications and another for incorrect classifications. It's a switch with two possible values: right, or wrong.
         Weighting scheme 3 is designed for projects where the gold-standards have variable intrinsic difficulty, so it sets the seed value for each gold-standard subject depending on how many people correctly classified it. Then it sums those individual seeds to the seed that's passed here and computes the weight in the same way as weighting scheme 2.
        '''

    elif (which_weight == 2) | (which_weight == 3):
        if n_gs < 1: # don't divide by or take the log of 0

            # also if they didn't do any gold-standard classifications assume they have the default weight
            return c0
        else:
            # note the max of 3 is unlikely to be reached, but someone could hit the floor.
            # can't get higher than 3 or lower than 0.05
            # return np.min([3.0, np.max([0.05, c0*pow((1.0 + np.log10(n_gs)), (float(seed)/float(n_gs)))])])
            # return np.min([3.0,(float(seed)/float(n_gs))])
            return np.min([3.0, np.max([0.05,(float(seed)/float(n_gs))])])
            
    else:
        # unweighted - so maybe don't even enter this function if which_weight is not 1 or 2...
        return 1.0


#################################################################################
#################################################################################
#################################################################################


# Get the Gini coefficient - https://en.wikipedia.org/wiki/Gini_coefficient
# Typical values of the Gini for healthy Zooniverse projects (Cox et al. 2015) are
#  in the range of 0.7-0.9.

def gini(list_of_values):
    sorted_list = sorted(list_of_values)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(list_of_values) / 2
    return (fair_area - area) / fair_area












