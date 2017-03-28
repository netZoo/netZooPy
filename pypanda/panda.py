from __future__ import print_function

import pandas as pd
import numpy as np
from scipy.stats import zscore

import time
import math
import sys

class Panda(object):
    '''Import and run PANDA algorithm.'''
    def __init__(self, expression_file, motif_file, ppi_file, gene_names=None):
        '''Load expression, motif and optional ppi data.'''
        # Create coexpression network
        t1 = time.time()
        print('Loading expression data ...')
        if expression_file.endswith('.npy'):
            self.expression_data = np.load(expression_file)
            self.gene_names = np.load(gene_names).tolist()
        elif expression_file.endswith('.pickle'):
            self.expression_data = pd.read_pickle(expression_file)
            self.gene_names = self.expression_data.index.tolist()
        else:
            self.expression_data = pd.read_table(expression_file, sep='\t', header=None, comment='#', index_col=0)
            self.gene_names = self.expression_data.index.tolist()
        self.num_genes = len(self.gene_names)
        print(self.expression_data.shape)

        # Create motif network
        print('Loading motif data ...')
        if motif_file.endswith('.npz'):
            self.motif_data = np.load(motif_file)
        elif motif_file.endswith('.pickle'):
            self.motif_data = pd.read_pickle(motif_file)
        else:
            self.motif_data = pd.read_table(motif_file, sep='\t', header=None, comment='#')
        self.unique_tfs = sorted(set(self.motif_data[0]))
        self.num_tfs = len(self.unique_tfs)
        print('Unique TFs:', self.num_tfs)

        print('Creating motif network ...')
        self.motif_matrix = np.zeros((self.num_tfs, self.num_genes))
        tf2idx = {x: i for i,x in enumerate(self.unique_tfs)}
        gene2idx = {x: i for i,x in enumerate(self.gene_names)}
        idx_tfs = [tf2idx[x] for x in self.motif_data[0]]
        idx_genes = [gene2idx[x] for x in self.motif_data[1]]
        idx = np.ravel_multi_index((idx_tfs, idx_genes), self.motif_matrix.shape)
        self.motif_matrix.ravel()[idx] = self.motif_data[2]
        #self.motif_matrix.ravel()[idx] = 1

        # Create PPI network
        print('Loading PPI data ...')
        if ppi_file.endswith('.npz'):
            self.ppi_data = np.load(ppi_file)
        elif ppi_file.endswith('.pickle'):
            self.ppi_data = pd.read_pickle(ppi_file)
        else:
            self.ppi_data = pd.read_table(ppi_file, sep='\t', header=None, comment='#')

        print('Creating PPI network ...')
        self.ppi_matrix = np.identity(self.num_tfs)
        idx_tf1 = [tf2idx[x] for x in self.ppi_data[0]]
        idx_tf2 = [tf2idx[x] for x in self.ppi_data[1]]
        idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
        self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
        idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
        self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
        print('Total elapsed time:', time.time() - t1, 'sec.')

        print('Calculating coexpression network ...')
        t2 = time.time()
        self.correlation_matrix = np.corrcoef(self.expression_data)
        print('Elapsed time:', time.time() - t2, 'sec.')

        print('Normalizing coexpression matrix ...')
        t3 = time.time()
        self.correlation_matrix = self.__normalize_network(self.correlation_matrix)
        print('Normalizing PPI network ...')
        self.ppi_matrix = self.__normalize_network(self.ppi_matrix)
        print('Normalizing motif network ...')
        self.motif_matrix = self.__normalize_network(self.motif_matrix, square=False)
        print('Total elapsed time:', time.time() - t3, 'sec.')

        #print('Saving coexpression network ...')
        #t1 = time.time()
        #np.save('/dev/shm/correlation_matrix.normalized.npy', self.correlation_matrix)
        #print('Elapsed time:', time.time() - t1, 'sec.')

        print('Saving motif network ...')
        t1 = time.time()
        np.save('/dev/shm/normalized_motif_matrix.npy', self.motif_matrix)
        #np.save('/dev/shm/unique_tfs.npy', self.unique_tfs)
        print('Elapsed time:', time.time() - t1, 'sec.')

        print('Saving PPI network ...')
        t1 = time.time()
        np.save('/dev/shm/normalized_ppi_matrix.npy', self.ppi_matrix)
        print('Elapsed time:', time.time() - t1, 'sec.')

        # Check if nan in the networks
        # if one column or row has all 0 values (or all the same values), stdev will be 0, leading to z = nan
        #assert not np.isnan(self.motif_matrix).any()
        #assert not np.isnan(self.ppi_matrix).any()
        #assert not np.isnan(self.correlation_matrix).any()

        # Run panda algorithm
        print('Running PANDA algorithm ...')
        self.panda_network = self.panda_loop(self.correlation_matrix, self.motif_matrix, self.ppi_matrix, step_print=True)

        #print('Saving PANDA network to /dev/shm ...')
        #t1 = time.time()
        #np.save('/dev/shm/panda_network.npy', self.panda_network)
        #print('Elapsed time:', time.time() - t1, 'sec.')

    def __normalize_network(self, x, square=True):
        if square:
            norm_col = zscore(x, axis=0)
            return (norm_col + norm_col.T) / math.sqrt(2)
        else:
            norm_col = zscore(x, axis=0)
            norm_row = zscore(x, axis=1)
            return (norm_col + norm_row) / math.sqrt(2)

    def panda_loop(self, correlation_matrix, motif_matrix, ppi_matrix, step_print=True):
        '''Run panda algorithm.'''
        def t_function(x, y=None):
            '''T function.'''
            if y is None:
                a_matrix = np.dot(x, x.T)
                s = np.square(x).sum(axis=1)
                a_matrix /= np.sqrt(np.tile(s, (x.shape[0],1)) + s.reshape(-1, 1) - np.abs(a_matrix))
            else:
                a_matrix = np.dot(x, y)
                a_matrix /= np.sqrt(np.tile(np.square(y).sum(axis=0), (x.shape[0],1)) + np.square(x).sum(axis=1).reshape(-1, 1) - np.abs(a_matrix))
            return a_matrix
        def update_diagonal(diagonal_matrix, num, alpha, step):
            '''Update diagonal.'''
            np.fill_diagonal(diagonal_matrix, np.nan)
            diagonal_std = np.nanstd(diagonal_matrix, 1)
            diagonal_fill = diagonal_std * num * math.exp(2 * alpha * step)
            np.fill_diagonal(diagonal_matrix, diagonal_fill)
        panda_loop_time = time.time()
        step = 0
        hamming = 1
        alpha = 0.1
        while hamming > 0.001:
            t = time.time()
            # Update motif_matrix
            W = 0.5 * (t_function(ppi_matrix, motif_matrix) + t_function(motif_matrix, correlation_matrix))  # W = (R + A) / 2
            hamming = np.abs(motif_matrix - W).mean()
            motif_matrix *= (1 - alpha)
            motif_matrix += (alpha * W)

            if hamming > 0.001:
                # Update ppi_matrix
                ppi = t_function(motif_matrix)  # t_func(X, X.T)
                update_diagonal(ppi, self.num_tfs, alpha, step)
                ppi_matrix *= (1 - alpha)
                ppi_matrix += (alpha * ppi)

                # Update correlation_matrix
                motif = t_function(motif_matrix.T)  # t_func(X.T, X)
                update_diagonal(motif, self.num_genes, alpha, step)
                correlation_matrix *= (1 - alpha)
                correlation_matrix += (alpha * motif)
            if step_print:
                print('step: {}, hamming: {}, elapsed: {:.2f} sec.'.format(step, hamming, time.time() - t))
            step = step + 1
        print('running panda took: %.2f seconds' % (time.time() - panda_loop_time))
        return motif_matrix

    def __panda_results_data_frame(self):
        '''Results to data frame.'''
        tfs = np.tile(self.unique_tfs, (len(self.gene_names), 1)).flatten()
        genes = np.tile(self.gene_names, (len(self.unique_tfs), 1)).transpose().flatten()
        motif = self.motif_matrix.transpose().flatten()
        force = self.panda_network.transpose().flatten()
        self.flat_panda_network = force
        self.export_panda_results = pd.DataFrame({'tf':tfs, 'gene': genes,'motif': motif, 'force': force})
        self.export_panda_results = self.export_panda_results[['tf', 'gene', 'motif', 'force']]

    def save_panda_results(self, file = 'panda.pairs'):
        '''Write results to file.'''
        self.__panda_results_data_frame()
        self.export_panda_results.to_csv(file, index=False, header=False, sep="\t")
