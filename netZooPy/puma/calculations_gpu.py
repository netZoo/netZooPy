from __future__ import print_function

import numpy as np
import cupy as cp

from netZooPy.panda.calculations_gpu import gt_function, gupdate_diagonal


def compute_puma_gpu(
    correlation_matrix,
    motif_matrix,
    ppi_matrix,
    sorted_index,
    alpha=0.1,
    threshold=0.001,
):

    """
    Description:
        The PUMA algorithm.

    Inputs:
        correlation_matrix: Input coexpression matrix.
        motif_matrix      : Input motif regulation prior network.
        ppi_matrix        : Input PPI matrix.
        computing         : 'cpu' uses Central Processing Unit (CPU) to run PANDA.
                            'gpu' use the Graphical Processing Unit (GPU) to run PANDA.

    Uses (from panda):
        t_function      : Continuous Tanimoto similarity function computed on the CPU.
        update_diagonal : Updates the diagonal of the input matrix in the message passing computed on the CPU.
        gt_function     : Continuous Tanimoto similarity function computed on the GPU.
        gupdate_diagonal: Updates the diagonal of the input matrix in the message passing computed on the GPU.
    """

    num_tfs, num_genes = motif_matrix.shape
    step = 0
    hamming = 1
    s1 = sorted_index.copy()

    # Alessandro
    TFCoopInit = ppi_matrix.copy()

    ppi_matrix = cp.array(ppi_matrix)
    motif_matrix = cp.array(motif_matrix)
    correlation_matrix = cp.array(correlation_matrix)

    while hamming > threshold:
        W = 0.5 * (
            gt_function(ppi_matrix, motif_matrix)
            + gt_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        hamming = cp.abs(motif_matrix - W).mean()
        motif_matrix = cp.array(motif_matrix)
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W

        ppi = gt_function(motif_matrix)  # t_func(X, X.T)
        ppi = gupdate_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi

        # Alessandro
        TFCoopDiag = ppi_matrix.diagonal()
        ppi_matrix[s1] = TFCoopInit[s1]
        ppi_matrix[:, s1] = TFCoopInit[:, s1]
        np.fill_diagonal(ppi_matrix, TFCoopDiag)

        # Update correlation_matrix
        motif = gt_function(motif_matrix.T)
        motif = gupdate_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif

        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1

    motif_matrix = cp.asnumpy(motif_matrix)

    return motif_matrix
