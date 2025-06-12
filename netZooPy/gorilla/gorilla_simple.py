from __future__ import print_function
import math
import time
from tkinter.messagebox import RETRY
import pandas as pd
from scipy.stats import zscore
import numpy as np
from netZooPy.panda.panda import Panda
import sys
import os
import pymc3 as pm


class Gorilla(Panda):
    """
    GenOmic Regulatory Interaction Learning with Latent Attributes
        1. Reading in input data (co-expression, motif prior table, TF PPI data)
        2. Specifying the full model with pymc3
        3. Estimating regulatory bipartite network using variational inference


    Parameters
    ----------

            coexpression_matrix : array
                gene coexpression matrix of size (g,g) where g=number of genes
            prior_matrix : array
                TF-gene regulatory network based on TF motifs and/or chromatin accessibility as a
                matrix of size (t,g), where g=number of genes and t=number of TFs
            ppi_matrix : array
                TF-TF protein interaction network as a matrix of size (t,t), where t=number of TFs
            sparsify: bool
                If True some output edges are pruned based on posterior predictive region
                else, estimated bipartite network is returned as is
            typeI_error: float
                typeI error rate used for pruning edges from the estimated network
            output_folder: str
                folder where to save the results

    References
    ----------
    .. [1]__

    Authors: Enakshi Saha
    """

    def __init__(
            self,
            coexpression_matrix,
            prior_matrix,
            ppi_matrix=None,
            sparsify=False,
            typeI_error=0.05,
            output_folder='./gorilla_simple/'
    ):
        """Intialize instance of Panda class and load data."""

        self.coexpression_matrix = coexpression_matrix
        self.prior_matrix = prior_matrix
        self.ppi_matrix = ppi_matrix
        self.sparsify=sparsify
        self.typeI_error=typeI_error
        self.output_folder = output_folder

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)


    ########################
    ### METHODS ############
    ########################

    def run_gorilla(self, precision='single'):

        """Ligress algorithm

        Args:
            keep_coexpression (bool, optional): whether to save coexpression network
            save_memory (bool, optional): whether to save the coexpression with the gene names
        """

        if precision == 'single':
            atype = 'float32'
        elif precision == 'double':
            atype = 'float64'
        else:
            sys.exit('Precision %s unknonwn' % str(precision))

        # let's sort the expression and ppi data
        self.expression_data = self.expression_data.loc[self.universe_genes, :].astype(atype)
        correlation_complete = self.expression_data.T.corr()
