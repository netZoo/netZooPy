from __future__ import print_function

import math
import numpy as np
from scipy.stats import zscore
import sys

#
# Calculation functions:
# These functions were defined in Panda, but are
# also shared by Puma and possibly others.
# We are going to start putting these here so they can be
# shared by all classes as they are independent from the class


def t_function(x, y=None):
    """
    Description:
        Continuous Tanimoto similarity function computed on the CPU.

    Inputs:
        x: First object to measure the distance from. If only this matrix is provided, then the distance is meausred between the columns of x.
        y: Second object to measure the distance to.

    Ouputs:
        a_matrix: Matrix containing the pairwise distances.
    """
    if y is None:
        a_matrix = np.dot(x, x.T)
        s = np.square(x).sum(axis=1)
        # here the all zero entries for the motifs become nan
        a_matrix /= np.sqrt(s + s.reshape(-1, 1) - np.abs(a_matrix))
    else:
        a_matrix = np.dot(x, y)
        a_matrix /= np.sqrt(
            np.square(y).sum(axis=0)
            + np.square(x).sum(axis=1).reshape(-1, 1)
            - np.abs(a_matrix)
        )
    return a_matrix


def update_diagonal(diagonal_matrix, num, alpha, step):
    """
    Description:
        Updates the diagonal of the input matrix in the message passing computed on the CPU
    Inputs:
        diagonal_matrix: Input diagonal matrix.
        num            : Number of rows/columns.
        alpha          : Learning rate.
        step           : The current step in the algorithm.
    """
    np.fill_diagonal(diagonal_matrix, np.nan)
    diagonal_std = np.nanstd(diagonal_matrix, axis=0, ddof=0)
    diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
    np.fill_diagonal(diagonal_matrix, diagonal_fill)
    return diagonal_matrix

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)

def compute_panda_cpu(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    threshold=0.001,
    alpha=0.1,
):
    """Panda network optimization

    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """
    motif_matrix = motif_matrix.copy()
    ppi_matrix = ppi_matrix.copy()
    correlation_matrix = correlation_matrix.copy()

    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1

    while hamming > threshold:
        W = 0.5 * (
            t_function(ppi_matrix, motif_matrix)
            + t_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        hamming = np.abs(motif_matrix - W).mean()
        # Update motif matrix
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W
        # Update ppi_matrix
        # When motif is all zero, this goes to nan
        ppi = t_function(motif_matrix)  # t_func(X, X.T)
        ppi = update_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi
        # Update correlation_matrix
        motif = t_function(motif_matrix.T)
        motif = update_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif
        # del W, ppi, motif  # release memory for next step
        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1

    if math.isnan(hamming):
        print('Warning: NaN value for Hamming distance')
    return motif_matrix


def compute_panda(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    computing="cpu",
    threshold=0.001,
    alpha=0.1,
):
    """Panda network optimization

    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        computing (str) : either cpu or gpu. Defaults to 'cpu'
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """

    if computing == "cpu":
        # Initialise W and hamming
        motif_matrix = compute_panda_cpu(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            alpha=alpha,
        )

    elif computing == "gpu":
        from netZooPy.panda.calculations_gpu import compute_panda_gpu

        motif_matrix = compute_panda_gpu(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            alpha=alpha,
        )
    else:
        sys.error("ERROR: %s is not an existing computing device" % str(computing))
    return motif_matrix


def normalize_network(x):
    """
    Description:
        Standardizes the input data matrices.

    Inputs:
        x     : Input adjacency matrix.

    Outputs:
        normalized_matrix: Standardized adjacency matrix.
    """
    norm_col = zscore(x, ddof=0, axis=0)
    if (x.shape[0] == x.shape[1]) and check_symmetric(x):
        norm_row = norm_col.T
    else:
        norm_row = zscore(x, ddof=0, axis=1)
    # Alessandro: replace nan values
    normalized_matrix = (norm_col + norm_row) / math.sqrt(2)
    norm_total = (x - np.mean(x)) / np.std(x, ddof=1)  # NB zscore(x) is not the same
    nan_col = np.isnan(norm_col)
    nan_row = np.isnan(norm_row)
    normalized_matrix[nan_col] = (norm_row[nan_col] + norm_total[nan_col]) / math.sqrt(
        2
    )
    normalized_matrix[nan_row] = (norm_col[nan_row] + norm_total[nan_row]) / math.sqrt(
        2
    )
    normalized_matrix[nan_col & nan_row] = (
        2 * norm_total[nan_col & nan_row] / math.sqrt(2)
    )
    return normalized_matrix
