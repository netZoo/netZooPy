from __future__ import print_function

import sys
import numpy as np
from netZooPy.panda.calculations import (
    t_function,
    update_diagonal
)


def compute_puma_cpu(
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

    ppi_matrix = ppi_matrix.copy()
    motif_matrix = motif_matrix.copy()
    correlation_matrix = correlation_matrix.copy()

    while hamming > threshold:
        W = 0.5 * (
                t_function(ppi_matrix, motif_matrix)
                + t_function(motif_matrix, correlation_matrix)
        )  # W = (R + A) / 2
        hamming = np.abs(motif_matrix - W).mean()
        motif_matrix *= 1 - alpha
        motif_matrix += alpha * W

        ppi = t_function(motif_matrix)  # t_func(X, X.T)
        ppi = update_diagonal(ppi, num_tfs, alpha, step)
        ppi_matrix *= 1 - alpha
        ppi_matrix += alpha * ppi

        # Alessandro
        TFCoopDiag = ppi_matrix.diagonal().copy()
        ppi_matrix[s1] = TFCoopInit[s1]
        ppi_matrix[:, s1] = TFCoopInit[:, s1]
        np.fill_diagonal(ppi_matrix, TFCoopDiag)

        # Update correlation_matrix
        motif = t_function(motif_matrix.T)  # t_func(X.T, X)
        motif = update_diagonal(motif, num_genes, alpha, step)
        correlation_matrix *= 1 - alpha
        correlation_matrix += alpha * motif

        print("step: {}, hamming: {}".format(step, hamming))
        step = step + 1

    return motif_matrix


def compute_puma(
    correlation_matrix,
    motif_matrix,
    ppi_matrix,
    sorted_index,
    computing="cpu",
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

    if computing == "cpu":
        motif_matrix = compute_puma_cpu(
            correlation_matrix,
            motif_matrix,
            ppi_matrix,
            sorted_index,
            alpha=alpha,
        )

    elif computing == "gpu":
        from netZooPy.puma.calculations_gpu import compute_puma_gpu
        motif_matrix = compute_puma_gpu(
            correlation_matrix,
            motif_matrix,
            ppi_matrix,
            sorted_index,
            alpha=alpha,
        )
    else:
        sys.error("ERROR: %s is not an existing computing device")

    return motif_matrix
