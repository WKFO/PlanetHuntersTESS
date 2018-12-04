import numpy as np
import pandas as pd


def aggregate_class(grp):

    # get all the information for each "group" - one image/subject (except not subject in this case)

    """
    This function does the aggregation - groups into yes/no.
    
    """

    # translate the group to a dataframe (some indexing etc is different)
    thegrp = pd.DataFrame(grp)

    # figure out what we're looping over below
    answers = thegrp.planet_classification.unique()  # 'No' and 'Yes'
    
    # aggregating is a matter of grouping by different answers and summing the counts/weights
    byans = thegrp.groupby('planet_classification')
    ans_ct_tot = byans['count'].aggregate('sum')
    ans_wt_tot = byans['weight'].aggregate('sum')


    # we want fractions eventually, so we need denominators
    count_tot    = np.sum(ans_ct_tot) # we could also do len(thegrp)
    weight_tot   = np.sum(ans_wt_tot)

    # okay, now we should have a series of counts for each answer, one for weighted counts, and
    # the total votes and weighted votes for this subject.
    # now loop through the possible answers and create the raw and weighted vote fractions
    # and save the counts as well.

    # this is a list for now and we'll make it into a series and order the columns later
    # all the things you want in the output file

    class_agg = {}
    class_agg['count_unweighted'] = count_tot
    class_agg['count_weighted']   = weight_tot
    class_agg['subject_type']     = thegrp.subject_type.unique()[0]
    class_agg['subject']          = thegrp.subject.unique()[0]
    class_agg['candidate']        = thegrp.candidate.unique()[0]
    class_agg['TIC_ID']           = thegrp.TIC_ID.unique()[0]
    class_agg['SNR']              = thegrp.SNR.unique()[0]
    class_agg['min_pix']          = thegrp.min_pix.unique()[0]
    class_agg['max_pix']          = thegrp.max_pix.unique()[0]
    class_agg['min_time']         = thegrp.min_time.unique()[0]
    class_agg['max_time']         = thegrp.max_time.unique()[0]


    for a in answers: # create a new entry in the dataframe of yes and no (weighted and unweighted)

        # ignore blank classifications
        # don't be that jerk who labels things with "p0" or otherwise useless internal indices.
        # Use the text of the response next to this answer choice in the project builder (but strip spaces)
 
        try:
            raw_frac_label  = ('p_'+a).replace(' ', '_')
            wt_frac_label   = ('p_'+a+'_weight').replace(' ', '_')

        except:
            # if we're here it's because all the classifications for this
            # subject are empty, i.e. they didn't answer the question, just
            # clicked done
            # in that case you're going to crash in a second, but at least
            # print out some information about the group before you do

            print("OOPS ")
            print(thegrp)
            print('-----------')
            print(a)
            print('+++++')
            print(ans_ct_tot,ans_wt_tot,count_tot,weight_tot)
            print('~~~~~~')
            print(answers)
            print(len(a))

        class_agg[raw_frac_label] = ans_ct_tot[a]/float(count_tot)
        class_agg[wt_frac_label]  = ans_wt_tot[a]/float(weight_tot)

    # oops, this is hard-coded so that there's Yes and No as answers - sorry to those trying to generalise
    col_order = ["TIC_ID", "candidate", "p_Yes", "p_No", "p_Yes_weight", "p_No_weight",
                 "count_unweighted", "count_weighted", "subject_type", "subject", "SNR", "min_pix", "max_pix", "min_time", "max_time"]



    return pd.Series(class_agg)[col_order]
