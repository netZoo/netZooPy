from __future__ import print_function

import math
import time
import pandas as pd
import numpy as np
from scipy.stats import zscore
from .timer import Timer

class Panda(object):
    """ Using PANDA to infer gene regulatory network.

    1. Reading in input data (expression data, motif prior, TF PPI data)
    2. Computing coexpression network
    3. Normalizing networks
    4. Running PANDA algorithm
    5. Writing out PANDA network

    Authors: cychen, davidvi
    """
    def __init__(self, expression_file, motif_file, ppi_file, save_tmp=True):
        # =====================================================================
        # Data loading
        # =====================================================================
        with Timer('Loading expression data ...'):
            self.expression_data = pd.read_table(expression_file, sep='\t', header=None, index_col=0)
            self.gene_names = self.expression_data.index.tolist()
            self.num_genes = len(self.gene_names)
            print('Expression matrix:', self.expression_data.shape)

        with Timer('Loading motif data ...'):
            self.motif_data = pd.read_table(motif_file, sep='\t', header=None)
            self.unique_tfs = sorted(set(self.motif_data[0]))
            self.num_tfs = len(self.unique_tfs)
            print('Unique TFs:', self.num_tfs)

        with Timer('Loading PPI data ...'):
            self.ppi_data = pd.read_table(ppi_file, sep='\t', header=None)
            print('Number of PPIs:', self.ppi_data.shape[0])

        # Auxiliary dicts
        gene2idx = {x: i for i,x in enumerate(self.gene_names)}
        tf2idx = {x: i for i,x in enumerate(self.unique_tfs)}

        # =====================================================================
        # Network construction
        # =====================================================================
        with Timer('Calculating coexpression network ...'):
            self.correlation_matrix = np.corrcoef(self.expression_data)
            if np.isnan(self.correlation_matrix).any():
                np.fill_diagonal(self.correlation_matrix, 1)
                self.correlation_matrix = np.nan_to_num(self.correlation_matrix)

        with Timer('Creating motif network ...'):
            self.motif_matrix = np.zeros((self.num_tfs, self.num_genes))
            idx_tfs = [tf2idx[x] for x in self.motif_data[0]]
            idx_genes = [gene2idx[x] for x in self.motif_data[1]]
            idx = np.ravel_multi_index((idx_tfs, idx_genes), self.motif_matrix.shape)
            self.motif_matrix.ravel()[idx] = self.motif_data[2]

        with Timer('Creating PPI network ...'):
            self.ppi_matrix = np.identity(self.num_tfs)
            idx_tf1 = [tf2idx[x] for x in self.ppi_data[0]]
            idx_tf2 = [tf2idx[x] for x in self.ppi_data[1]]
            idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
            self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
            idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
            self.ppi_matrix.ravel()[idx] = self.ppi_data[2]

        # Clean up useless variables to release memory
        del self.motif_data, self.ppi_data, self.unique_tfs, self.gene_names

        # =====================================================================
        # Network normalization
        # =====================================================================
        with Timer('Normalizing networks ...'):
            self.correlation_matrix = self._normalize_network(self.correlation_matrix)
            self.motif_matrix = self._normalize_network(self.motif_matrix, square=False)
            self.ppi_matrix = self._normalize_network(self.ppi_matrix)

        # =====================================================================
        # Saving middle data to tmp
        # =====================================================================
        if save_tmp:
            with Timer('Saving expression matrix and normalized networks ...'):
                np.save('/tmp/expression.npy', self.expression_data.values)
                np.save('/tmp/motif.normalized.npy', self.motif_matrix)
                np.save('/tmp/ppi.normalized.npy', self.ppi_matrix)

        # Clean up useless variables to release memory
        del self.expression_data

        # =====================================================================
        # Running PANDA algorithm
        # =====================================================================
        print('Running PANDA algorithm ...')
        self.panda_network = self.panda_loop(self.correlation_matrix, self.motif_matrix, self.ppi_matrix)


    def _normalize_network(self, x, square=True):
        if square:
            norm_col = zscore(x, axis=0)
            return (norm_col + norm_col.T) / math.sqrt(2)
        else:
            norm_col = zscore(x, axis=0)
            norm_row = zscore(x, axis=1)
            return (norm_col + norm_row) / math.sqrt(2)

    def panda_loop(self, correlation_matrix, motif_matrix, ppi_matrix):
        """Panda algorithm.
        """
        def t_function(x, y=None):
            '''T function.'''
            if y is None:
                a_matrix = np.dot(x, x.T)
                s = np.square(x).sum(axis=1)
                a_matrix /= np.sqrt(s + s.reshape(-1, 1) - np.abs(a_matrix))
            else:
                a_matrix = np.dot(x, y)
                a_matrix /= np.sqrt(np.square(y).sum(axis=0) + np.square(x).sum(axis=1).reshape(-1, 1) - np.abs(a_matrix))
            return a_matrix

        def update_diagonal(diagonal_matrix, num, alpha, step):
            '''Update diagonal.'''
            np.fill_diagonal(diagonal_matrix, np.nan)
            diagonal_std = np.nanstd(diagonal_matrix, 1)
            diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
            np.fill_diagonal(diagonal_matrix, diagonal_fill)

        panda_loop_time = time.time()
        num_tfs, num_genes = motif_matrix.shape
        step = 0
        hamming = 1
        alpha = 0.1
        while hamming > 0.001:
            # Update motif_matrix
            W = 0.5 * (t_function(ppi_matrix, motif_matrix) + t_function(motif_matrix, correlation_matrix))  # W = (R + A) / 2
            hamming = np.abs(motif_matrix - W).mean()
            motif_matrix *= (1 - alpha)
            motif_matrix += (alpha * W)

            if hamming > 0.001:
                # Update ppi_matrix
                ppi = t_function(motif_matrix)  # t_func(X, X.T)
                update_diagonal(ppi, num_tfs, alpha, step)
                ppi_matrix *= (1 - alpha)
                ppi_matrix += (alpha * ppi)

                # Update correlation_matrix
                motif = t_function(motif_matrix.T)  # t_func(X.T, X)
                update_diagonal(motif, num_genes, alpha, step)
                correlation_matrix *= (1 - alpha)
                correlation_matrix += (alpha * motif)

                del W, ppi, motif  # release memory for next step

            print('step: {}, hamming: {}'.format(step, hamming))
            step = step + 1

        print('Running panda took: %.2f seconds!' % (time.time() - panda_loop_time))
        return motif_matrix

    def save_panda_results(self, path='panda.npy'):
        with Timer('Saving PANDA network to %s ...' % path):
            if path.endswith('.txt'):
                np.savetxt(path, self.panda_network)
            elif path.endswith('.csv'):
                np.savetxt(path, self.panda_network, delimiter=',')
            elif path.endswith('.tsv'):
                np.savetxt(path, self.panda_network, delimiter='/t')
            else:
                np.save(path, self.panda_network)
