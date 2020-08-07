from __future__ import print_function

import os, os.path,sys
import numpy as np
sys.path.insert(1,'../panda')
from netZooPy.puma.puma import Puma
from .timer import Timer

class LionessPuma(Puma):
    """
    Description:
         Using LIONESS to infer single-sample gene regulatory networks.

    Usage:
        1. Reading in PUMA network and preprocessed middle data
        2. Computing coexpression network
        3. Normalizing coexpression network
        4. Running PUMA algorithm
        5. Writing out LIONESS networks

    Authors: 
        cychen, davidvi
    """

    def __init__(self, obj, start=1, end=None, save_dir='lioness_output', save_fmt='npy'):
        # Load data
        with Timer("Loading input data ..."):
            self.expression_matrix = obj.expression_matrix
            self.motif_matrix = obj.motif_matrix
            self.ppi_matrix = obj.ppi_matrix
            
            tfs = np.tile(obj.unique_tfs, (len(obj.gene_names), 1)).flatten()
            genes = np.repeat(obj.gene_names, obj.num_tfs)
            motif = obj.motif_matrix_unnormalized.flatten(order='F')
            self.export_lioness_results = np.column_stack((tfs,genes,motif))
            
            if hasattr(obj,'puma_network'):
                self.network = obj.puma_network
            elif hasattr(obj,'puma_network'):
                self.network = obj.puma_network
            else:
                print('Cannot find puma network in object')
                raise AttributeError('Cannot find puma network in object')
            # Alessandro
            self.s1 = obj.s1
            del obj

        # Get sample range to iterate
        self.n_conditions = self.expression_matrix.shape[1]
        self.indexes = range(self.n_conditions)[start-1:end]  # sample indexes to include

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Run LIONESS
        self.__lioness_loop()

        # create result data frame
        #self.export_lioness_results = pd.DataFrame(self.lioness_network)

    def __lioness_loop(self):
        for i in self.indexes:
            print("Running LIONESS for sample %d:" % (i+1))
            idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i

            with Timer("Computing coexpression network:"):
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

            with Timer("Normalizing networks:"):
                correlation_matrix = self._normalize_network(correlation_matrix)

            with Timer("Inferring LIONESS network:"):
                subset_puma_network = self.puma_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix))
                lioness_network = self.n_conditions * (self.network - subset_puma_network) + subset_puma_network
            
            force = lioness_network.flatten(order='F')
            self.export_lioness_results = np.column_stack((self.export_lioness_results, force))
                
        return 

    def save_lioness_results(self, path='lioness.txt'):
        '''Write lioness results to file.'''
        #self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        if path.endswith('.txt'):
            np.savetxt(path, self.export_lioness_results, fmt='%s', delimiter=" ", header="")
        elif path.endswith('.csv'):
            np.savetxt(path, self.export_lioness_results, fmt='%s', delimiter=",", header="")
        elif path.endswith('.tsv'):
            np.savetxt(path, self.export_lioness_results, fmt='%s', delimiter="\t", header="")
        else:
            np.save(path, self.export_lioness_results)
        return None
