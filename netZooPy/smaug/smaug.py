from __future__ import print_function
import math
import time
import pandas as pd
from netZooPy.smaug.timer import Timer
import numpy as np
import sys
import os
import pandas as pd
import scipy.stats as stats
from pylab import meshgrid
from scipy.optimize import minimize
from scipy import optimize
import scipy.special as sc
import scipy.integrate as integrate
import statsmodels.stats.multitest as multi
from netZooPy.dragon import estimate_penalty_parameters_dragon
from netZooPy.dragon import get_shrunken_covariance_dragon


def compute_smaug(expression_matrix, methylation_matrix, expression_mean, methylation_mean, sample_idx, online_coexpression=False, computing='cpu', cores=1,
                   delta=None, compute_sparse=False, confidence=0.05, save_pvals=False):
    """Compute one smaug matrix. Takes as input an expression matrix, a methylation matrix, and the
    index of the sample for which smaug is computed (index of the columns).

    expression_matrix (numpy array): genes as rows, samples as columns
    expression_mean (numpy array): row-wise mean of expression.
    methylation_matrix (numpy array): methylation probes as rows, samples as columns
    methylation_mean (numpy array): row-wise mean of methylation.
    sample_idx (int): index of the sample of interest
    delta (float, optional): delta value for the computation. Defaults to None, which means that delta is tuned

    """

    mask_include = [True] * expression_matrix.shape[1]
    mask_include[sample_idx] = False

    # Compute covariance matrix from the rest of the data, leaving out sample
    lambdas = estimate_penalty_parameters_dragon(expression_matrix[:, mask_include], methylation_matrix[:, mask_include])
    covariance_matrix = get_shrunken_covariance_dragon(expression_matrix[:, mask_include], methylation_matrix[:, mask_include], lambdas)

    # Compute posterior weight delta from data
    if (delta == None):
        delta = 1 / (3 + 2 * np.sqrt(covariance_matrix.diagonal()).mean() / covariance_matrix.diagonal().var())
    else:
        assert type(delta) == float

    # Append expression amd methylation matrix by row
    combined_data = np.concatenate((expression_matrix, methylation_matrix), axis=0)

    # Append mean expression and mean methylation
    combined_mean = np.concatenate((expression_mean, methylation_mean))

    # Compute sample-specific covariance matrix
    sscov = delta * np.outer((combined_data - combined_mean)[:, sample_idx],
                             (combined_data - combined_mean)[:, sample_idx]) + (1 - delta) * covariance_matrix

    # Compute sample-specific DRAGON from sample-specific covariance
    Theta = np.linalg.inv(sscov)
    p = Theta.shape[0]
    A = np.sqrt(np.zeros((p, p)) + np.diag(Theta))
    smaug_matrix = -Theta / A / A.T
    smaug_matrix = smaug_matrix - np.diag(np.diag(smaug_matrix))

    return (smaug_matrix)





