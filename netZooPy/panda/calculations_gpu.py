from __future__ import print_function

import math
import cupy as cp
import numpy as np
#below are the imports for memory and time profiling
#import time
#remove this
# importing the library
#from memory_profiler import profile

#
# GPU Calculation functions:
# These functions were defined in Panda, but are
# also shared by Puma and possibly others.
# We are going to start putting these here so they can be
# shared by all classes as they are independent from the class

def gt_function(x, y=None):
    """
    Description:
        Continuous Tanimoto similarity function computed on the GPU
    Inputs:
        x: First object to measure the distance from. If only this matrix is provided, then the distance is meausred between the columns of x.
        y: Second object to measure the distance to
    Ouputs:
        a_matrix: Matrix containing the pairwsie distances.
    """
    if y is None:
        a_matrix = cp.dot(x, x.T)
        s = cp.square(x).sum(axis=1)
        a_matrix /= cp.sqrt(s + s.reshape(-1, 1) - cp.abs(a_matrix))
        #del x
    else:
        a_matrix = cp.dot(x, y)
        a_matrix /= cp.sqrt(
            cp.square(y).sum(axis=0)
            + cp.square(x).sum(axis=1).reshape(-1, 1)
            - cp.abs(a_matrix)
        )
        #del x,y
        
    #cp._default_memory_pool.free_all_blocks()
    return a_matrix

def gupdate_diagonal(diagonal_matrix, num, alpha, step):
    """
    Description:
        Updates the diagonal of the input matrix in the message passing computed on the GPU
    Inputs:
        diagonal_matrix: Input diagonal matrix.
        num            : Number of rows/columns.
        alpha          : Learning rate.
        step           : The current step in the algorithm.
    """
    cp.fill_diagonal(diagonal_matrix, cp.nan)
    #diagonal_std = cp.nanstd(diagonal_matrix, 1)
    diagonal_std = cp.nanstd(diagonal_matrix, axis=0, ddof=0)
    diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
    cp.fill_diagonal(diagonal_matrix, diagonal_fill)
    #del diagonal_fill
    #cp._default_memory_pool.free_all_blocks()
    
    return diagonal_matrix

# Use this for memory profiling
#@profile
def compute_panda_gpu(
    correlation_matrix,
    ppi_matrix,
    motif_matrix,
    threshold=0.001,
    alpha=0.1,
):

    """
    Panda network optimization

    Args:
        correlation_matrix (numpy float): coexpression matrix
        ppi_matrix (numpy float): PPI network matrix
        motif_matrix (numpy float): motif matrix
        threshold (float, optional): hamming distance threshold for stop. Defaults to 0.001.
        alpha (float, optional): learning rate. Defaults to 0.1
    """
    print("Computing panda on GPU")
    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1

    ppi_matrix = cp.array(ppi_matrix.copy())
    motif_matrix = cp.array(motif_matrix.copy())
    correlation_matrix = cp.array(correlation_matrix.copy())
    
    while hamming > threshold:

        W = 0.5 * (
            gt_function(ppi_matrix, motif_matrix)
            + gt_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        #print(W.dtype)
        #print(ppi_matrix.dtype)
        #print(motif_matrix.dtype)
        hamming = cp.abs(motif_matrix - W).mean()
        #print(hamming.dtype)
        # update motif matrix
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W
        # Update ppi_matrix
        ppi = gt_function(motif_matrix)  # t_func(X, X.T)
        ppi = gupdate_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi

        # Update correlation_matrix
        motif = gt_function(motif_matrix.T)
        motif = gupdate_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif
        # del W, ppi, motif  # release memory for next step
        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1
        
        del W, ppi, motif  # release memory for next step

    if math.isnan(hamming):
        print('Warning: NaN value for Hamming distance')

    motif_matrix_np = cp.asnumpy(motif_matrix)
    del motif_matrix
    cp._default_memory_pool.free_all_blocks()
    
    return motif_matrix_np
