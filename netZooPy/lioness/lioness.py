from __future__ import print_function

import os, os.path,sys
import numpy as np
import pandas as pd
from .timer import Timer
sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda

class Lioness(Panda):
    """
    Description:
       Using LIONESS to infer single-sample gene regulatory networks.

    Usage:
       1. Reading in PANDA network and preprocessed middle data
       2. Computing coexpression network
       3. Normalizing coexpression network
       4. Running PANDA algorithm
       5. Writing out LIONESS networks

    Inputs:
       obj: PANDA object, generated with keep_expression_matrix=True.
       obj.motif_matrix: TF DNA motif binding data in tf-by-gene format.
                         If set to None, Lioness will be performed on gene coexpression network.
       computing  : 'cpu' uses Central Processing Unit (CPU) to run PANDA
                    'gpu' use the Graphical Processing Unit (GPU) to run PANDA

    Authors: 
       cychen, davidvi
    """

    def __init__(self, obj, computing='cpu', precision='double',start=1, end=None, save_dir='lioness_output', save_fmt='npy'):
        # Load data
        with Timer("Loading input data ..."):
            self.export_panda_results = obj.export_panda_results
            self.expression_matrix = obj.expression_matrix
            self.motif_matrix = obj.motif_matrix
            self.ppi_matrix = obj.ppi_matrix
            self.correlation_matrix=obj.correlation_matrix
            if precision=='single':
                self.correlation_matrix=np.float32(self.correlation_matrix)
                self.motif_matrix=np.float32(self.motif_matrix)
                self.ppi_matrix=np.float32(self.ppi_matrix)
            self.computing=computing
            if hasattr(obj,'panda_network'):
                self.network = obj.panda_network
            elif hasattr(obj,'puma_network'):
                self.network = obj.puma_network
            else:
                print('Cannot find panda or puma network in object')
                raise AttributeError('Cannot find panda or puma network in object')
            del obj

        # Get sample range to iterate
        self.n_conditions = self.expression_matrix.shape[1]
        self.indexes = range(self.n_conditions)[start-1:end]  # sample indexes to include

        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        # prepare for onlineCoexpression
        self.fullmean=np.mean(self.expression_matrix,axis=1,dtype=np.float64)
        self.fullstd=np.std(self.expression_matrix,axis=1,dtype=np.float64)
        self.fullcov=np.cov(self.expression_matrix)
        # Run LIONESS
        self.total_lioness_network = self.__lioness_loop()

        # create result data frame
        self.export_lioness_results = pd.DataFrame(self.total_lioness_network)

    def __lioness_loop(self):
        for i in self.indexes:
            print("Running LIONESS for sample %d:" % (i+1))
            idx = [x for x in range(self.n_conditions) if x != i]  # all samples except i
            with Timer("Computing coexpression network:"):
                subj_exp=self.expression_matrix[:, i]
                if self.computing=='gpu':
                    import cupy as cp
                    correlation_matrix = cp.corrcoef(self.expression_matrix[:, idx])
                    if cp.isnan(correlation_matrix).any():
                        cp.fill_diagonal(correlation_matrix, 1)
                        correlation_matrix = cp.nan_to_num(correlation_matrix)
                    correlation_matrix=cp.asnumpy(correlation_matrix)
                else:
                    correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                    if np.isnan(correlation_matrix).any():
                        np.fill_diagonal(correlation_matrix, 1)
                        correlation_matrix = np.nan_to_num(correlation_matrix)

            with Timer("Normalizing networks:"):
                correlation_matrix_orig = correlation_matrix # save matrix before normalization
                correlation_matrix = self._normalize_network(correlation_matrix)

            with Timer("Inferring LIONESS network:"):
                if self.motif_matrix is not None:
                    del correlation_matrix_orig
                    subset_panda_network = self.panda_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix),self.computing)
                else:
                    del correlation_matrix
                    subset_panda_network = correlation_matrix_orig

            lioness_network = self.n_conditions * (self.network - subset_panda_network) + subset_panda_network

            with Timer("Saving LIONESS network %d to %s using %s format:" % (i+1, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "lioness.%d.%s" % (i+1, self.save_fmt))
                if self.save_fmt == 'txt':
                    np.savetxt(path, lioness_network)
                elif self.save_fmt == 'npy':
                    np.save(path, lioness_network)
                elif self.save_fmt == 'mat':
                    from scipy.io import savemat
                    savemat(path, {'PredNet': lioness_network})
                else:
                    print("Unknown format %s! Use npy format instead." % self.save_fmt)
                    np.save(path, lioness_network)
            if i == 0:
                self.total_lioness_network = np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)
            else:
                self.total_lioness_network=np.column_stack((self.total_lioness_network ,np.fromstring(np.transpose(lioness_network).tostring(),dtype=lioness_network.dtype)))

        return self.total_lioness_network

    def save_lioness_results(self, file='lioness.txt'):
        '''Write lioness results to file.'''
        #self.lioness_network.to_csv(file, index=False, header=False, sep="\t")
        np.savetxt(file, self.total_lioness_network, delimiter="\t",header="")
        return None

