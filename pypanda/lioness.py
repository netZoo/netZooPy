from __future__ import print_function

from .panda import Panda

import time
import pandas as pd
import numpy as np
import functools
import math

class Lioness(Panda):
    def __init__(self, panda_data):
        '''Import values from panda for lioness.'''
        self.export_panda_results = panda_data.export_panda_results
        self.expression_matrix = panda_data.expression_matrix
        self.motif_data = panda_data.motif_data
        if self.motif_data is not None:
            self.motif_matrix = panda_data.motif_matrix
            self.ppi_matrix = panda_data.ppi_matrix
            self.num_genes = panda_data.num_genes
            self.num_tfs = panda_data.num_tfs
        self.panda_network = panda_data.panda_network
        self.expression_data = panda_data.expression_data
        self.flat_panda_network = panda_data.flat_panda_network
        #run lioness
        self.__lioness()
        #create result data frame
        self.__lioness_results_data_frame()
        return None
    def __lioness(self):
        '''Run lioness on network.'''
        def lioness_loop(condition, number_conditions, flat_panda_network, expression_matrix, motif_matrix, ppi_matrix):
            '''Lioness algorithm.'''
            idx = range(0, number_conditions)
            idx.remove(condition)
            subset_expression_matrix = expression_matrix[:,idx]
            correlation_matrix = np.corrcoef(subset_expression_matrix)
            if self.motif_data is not None:
                subset_panda_network = self.panda_loop(correlation_matrix, motif_matrix, ppi_matrix)
            else:
                subset_panda_network = correlation_matrix
            subset_panda_network = subset_panda_network.transpose().flatten()
            lioness_network = number_conditions*(flat_panda_network-subset_panda_network)+subset_panda_network
            return lioness_network
        lioness_loop_time = time.time()
        number_conditions = self.expression_matrix.shape[1]
        if self.motif_data is not None:
            self.lioness_network = np.zeros((self.num_genes*self.num_tfs, number_conditions))
        else:
            self.lioness_network = np.zeros((self.expression_matrix.shape[0]*self.expression_matrix.shape[0], number_conditions))
        if number_conditions < 4:
            print('Running lioness requires at least 4 samples.')
            return None
        for condition in range(0, number_conditions):
            if self.motif_data is not None:
                self.lioness_network[:,condition] = lioness_loop(condition = condition,
                                                            number_conditions = number_conditions,
                                                            flat_panda_network = self.flat_panda_network,
                                                            expression_matrix = self.expression_matrix,
                                                            motif_matrix = self.motif_matrix,
                                                            ppi_matrix = self.ppi_matrix)
            else:
                self.lioness_network[:,condition] = lioness_loop(condition = condition,
                                                            number_conditions = number_conditions,
                                                            flat_panda_network = self.flat_panda_network,
                                                            expression_matrix = self.expression_matrix,
                                                            motif_matrix = None,
                                                            ppi_matrix = None)
        self.lioness_network = np.matrix(self.lioness_network)
        print('running lioness took: %s seconds' % (time.time() - lioness_loop_time))
        return None
    def __lioness_results_data_frame(self):
        ''''Results to data frame.'''
        self.export_lioness_results = pd.DataFrame(self.lioness_network)
        return None
    def save_lioness_results(self, file = 'lioness.txt'):
        '''Write lioness results to file.'''
        self.export_lioness_results.to_csv(file, index=False, header=False, sep="\t")
        return None
