from __future__ import print_function

import os, os.path
import numpy as np
from .panda import Panda
from .timer import Timer

class Lioness(Panda):
    """Using LIONESS to infer single-sample gene regulatory networks.

    1. Reading in PANDA network and preprocessed middle data
    2. Computing coexpression network
    3. Normalizing coexpression network
    4. Running PANDA algorithm
    5. Writing out LIONESS networks

    Authors: cychen, davidvi
    """
    def __init__(self, expression_matrix, motif_matrix, ppi_matrix, panda_network, start=1, end=None, save_dir='lioness_output', save_fmt='npy'):
        # Load data
        with Timer("Loading input data ..."):
            self.expression_matrix = np.load(expression_matrix)
            self.motif_matrix = np.load(motif_matrix)
            self.ppi_matrix = np.load(ppi_matrix)
            self.panda_network = np.load(panda_network)

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

    def __lioness_loop(self):
        for i in self.indexes:
            print("Running LIONESS for sample %d:" % (i+1))
            idx = [x for x in xrange(self.n_conditions) if x != i]  # all samples except i

            with Timer("Computing coexpression network:"):
                correlation_matrix = np.corrcoef(self.expression_matrix[:, idx])
                if np.isnan(correlation_matrix).any():
                    np.fill_diagonal(correlation_matrix, 1)
                    correlation_matrix = np.nan_to_num(correlation_matrix)

            with Timer("Normalizing networks:"):
                correlation_matrix = self._normalize_network(correlation_matrix)

            with Timer("Inferring LIONESS network:"):
                subset_panda_network = self.panda_loop(correlation_matrix, np.copy(self.motif_matrix), np.copy(self.ppi_matrix))
                lioness_network = self.n_conditions * (self.panda_network - subset_panda_network) + subset_panda_network

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
