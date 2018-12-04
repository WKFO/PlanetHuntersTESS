import numpy as np
import pandas as pc

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from sklearn import metrics
from sklearn.metrics import roc_auc_score

outpath = '/Users/Nora/Documents/research/TESS/planethunters/Rel01/output/figures/'


def ROC_curve_new(data, weight = None):


	true = data['subject_type'] # array of true and false where true has a planet or EB

	true_scores  = data['p_Yes']
	true_scores_w  = data['p_Yes_weight']

	fpr, tpr, thresholds = metrics.roc_curve(true, true_scores, pos_label=True, drop_intermediate = False)
	fpr_w, tpr_w, thresholds_w = metrics.roc_curve(true, true_scores_w, pos_label=True, drop_intermediate = False)

	print(fpr, tpr, thresholds)

	#print (len(thresholds))
	fig,axis = plt.subplots(figsize=(9,6))
	
	plt.plot(fpr, tpr, lw=3, alpha=0.9)
	plt.plot(fpr_w, tpr_w, lw=3, alpha=0.9)

	plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
	
	plt.xlim([-0.05, 1.05])
	plt.ylim([-0.05, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")

	print ("saving figures to file...")
	plt.savefig(outpath + '/ROC_curve_sims.png', format='png', dpi=80,bbox_inches='tight')






