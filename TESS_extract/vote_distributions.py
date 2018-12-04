
import numpy as np
from scipy.stats import norm


def get_net_voting_distribution(xvals, widths, weight, x_eval):
    """
    Find smooth pdf of contribution of user

    Args:
        xvals (list): x's of transits identified by user
        widths (list): widths of transits identified by user
        weight (float): relative weight of user, for Gaussian area

    Returns:
        (np.array): net voting distribution, evaluated at x_eval points
    """
    assert len(xvals) == len(widths)

    classification_pdf = np.zeros_like(x_eval)
    for classification_n in range(len(xvals)):
        x = xvals[classification_n]
        width = widths[classification_n]
        classification_pdf += norm.pdf(x_eval, loc=x, scale=width/2.355)   # divide the width by two to use the FWHM instead of the standdard deviation

    return classification_pdf * weight  # scale up by user
