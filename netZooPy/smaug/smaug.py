from __future__ import print_function
import time
from netZooPy.smaug.timer import Timer
import sys
import os
import pandas as pd
import numpy as np
import math
from netZooPy.smaug import io
from netZooPy.dragon import *      # To load DRAGON
from scipy.stats import norm # To get normal quantiles
from scipy.stats import false_discovery_control as fdr # To get adjusted p-values

class Smaug():
    """
    SMAUG


    Parameters
    ----------

            expression_file : str
                Path to file containing the gene expression data.
            methylation_file : str
                Path to file containing the methylation data.
            output_folder: str
                folder where to save the results
            delta: float
                posterior weight between 0 and 1 (If None (default) delta is tuned empirically from data)

    Notes
    ------

    Toy data:The example gene expression and methylation data that we have available here contains
    gene expression and methylation profiles for different samples in the columns.
    This is a small simulated excample.
    We provided these "toy" data so that the user can test the method.


    Sample SMAUG results:\b
        - Node1   Node2   Weight\n
        - gene1 cpg1	0.0	-0.951416589143\n
        - gene1 cpg2	0.0	-0.904241609324\n
        .
        .
        .
        - gene2 cpg1	0.0	-0.951416589143\n
        - gene2 cpg2	0.0	-0.904241609324\n

    Authors: Enakshi Saha
    """

    def __init__(
            self,
            expression_file,
            methylation_file
    ):
        """Intialize instance of Smaug class and load data."""

        self.expression_file = expression_file
        self.methylation_file = methylation_file

        # data read
        self.samples = None
        self.n_samples = None
        self.expression_data = None
        self.expression_genes = None
        self.expression_samples = None
        self.methylation_data = None
        self.methylation_probes = None
        self.methylation_samples = None
        # prepare all the data
        print('SMAUG: preparing expression and methylation')
        self._prepare_data()
        self.delta = None
        self.smaugs = []
        self.precisions = []
        self.pvals = []
        self.adjPvals = []

    ########################
    ### METHODS ############
    ########################
    def _prepare_data(self):
        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_data(
                self.expression_file, samples=self.samples
            )

        with Timer("Reading methylation data..."):
            # Read expression
            self.methylation_data, self.methylation_probes = io.prepare_data(
                self.methylation_file, samples=self.samples
            )

            self.expression_samples = self.expression_data.columns.tolist()
            self.methylation_samples = self.methylation_data.columns.tolist()

            self.expression_data = self.expression_data.T
            self.methylation_data = self.methylation_data.T

    def run_smaug(self, keep_in_memory=False, output_fmt=".hdf", output_folder='./smaug_output/',
                   delta=None, precision='single', sample_names=[]):
        """SMAUG algorithm

        Args:
            output_folder (str, optional): output folder. If an empty string is passed the matrix is automatically kept
            in memory, overwriting the value of keep_in_memory
            output_fmt (str, optional): format of output matrix. By default it is an hdf file, can be a txt or csv.
            keep_in_memory (bool, optional): if True, the partial correlation matrix is kept in memory, otherwise it is
            discarded after saving.
            only when the number of genes is not very big to avoid saving huge matrices.
            delta (float, optional): delta parameter. If default (None) delta is trained, otherwise pass a value.
            precision (str, optional): matrix precision, defaults to single precision.
        """

        smaug_start = time.time()

        # first let's reorder the expression data

        if precision == 'single':
            atype = 'float32'
        elif precision == 'double':
            atype = 'float64'
        else:
            sys.exit('Precision %s unknonwn' % str(precision))

        # sort expression and methylation data
        self.expression_data = self.expression_data.astype(atype)
        self.methylation_data = self.methylation_data.astype(atype)
        self.output_fmt = output_fmt
        self.output_folder = output_folder
        # If output folder is an empty string, keep the matrix in memory and don't save it to disk
        # Otherwise the output folder can be created and the matrix saved
        if self.output_folder == '':
            keep_in_memory = True
        else:
            if not os.path.exists(self.output_folder):
                os.makedirs(self.output_folder)

        # Compute lambda penalty parameters for DRAGON
        lambdas, lambdas_landscape = estimate_penalty_parameters_dragon(self.expression_data, self.methylation_data)

        # Compute population dragon
        pop_precision, _ = get_precision_matrix_dragon(self.expression_data, self.methylation_data, lambdas)

        # Append both omics data
        fulldata = np.append(self.expression_data, self.methylation_data, axis=1).T

        # Compute centered omics data
        z = fulldata - np.mean(fulldata, axis=0)

        print('SMAUG: We are starting to compute the networks...')
        if sample_names == []:
            sample_names = self.expression_samples
            sample_names = set(sample_names).intersection(set(self.methylation_samples))
        else:
            different = set(sample_names).difference(set(self.expression_samples).union(set(self.methylation_samples)))
            sample_names = set(sample_names).intersection(set(self.expression_samples))
            sample_names = set(sample_names).intersection(set(self.methylation_samples))
            if len(different) > 0:
                print('WARNING: some of the sample names are not in the expression data')
                print('\tMissing:')
                print('\t' + str(different))
                print('\tUsing:')
                print('\t' + str(sample_names))

        for s, sample in enumerate(sample_names):
            sample_start = time.time()
            # first run Smaug
            print('SMAUG: network for sample %s' % str(sample))
            if keep_in_memory:
                result_smaug, result_precision, result_pval_precision, result_pval_precision_adjusted = self.compute_individual_smaug(
                    fulldata, pop_precision, z, s, sample, delta, keep_in_memory)
                self.smaugs.append(result_smaug)
                self.precisions.append(result_precision)
                self.pvals.append(result_pval_precision)
                self.adjPvals.append(result_pval_precision_adjusted)
            else:
                self.compute_individual_smaug(fulldata, pop_precision, z, s, sample, delta, keep_in_memory)

        if keep_in_memory:
            return self

    def compute_individual_smaug(self, fulldata, pop_precision, z, s, sample, delta, keep_in_memory):
        """Runs smaug on one sample. All samples are saved separately.

        Args:
            fulldata: combined expression and methylation
            pop_precision: population level precision matrix
            z: centered omics data for the sample
            s, sample: sample index and name
            delta: delta parameter
            output_folder (str, optional): _description_.
        """

        mask_include = [True] * fulldata.shape[1]
        mask_include[s] = False

        print('SMAUG: computing network for sample %s' % str(sample))
        # Compute covariance matrix from the rest of the data, leaving out sample
        covariance_matrix = np.cov(fulldata[:, mask_include])

        # Compute posterior weight delta from data
        if delta == None:
            delta = 1 / (3 + 2 * np.sqrt(covariance_matrix.diagonal()).mean()
                         / covariance_matrix.diagonal().var())
        else:
            assert type(delta) == float

        # Compute sample-specific precision matrix
        M = np.outer(z[:, s], z[:, s])
        df = 1 / delta + covariance_matrix.shape[0] + 1
        numerator = pop_precision @ M @ pop_precision * delta / (1 - delta) ** 2
        denominator = 1 + delta / (1 - delta) * np.inner(z[:, s], pop_precision @ z[:, s])
        ssprecision = pop_precision / (1 - delta) - numerator / denominator

        # Compute sample-specific partial correlation
        p = ssprecision.shape[0]
        A = np.sqrt(np.zeros((p, p)) + np.diag(ssprecision))
        ssdragon = -ssprecision / A / A.T
        ssdragon = ssdragon - np.diag(np.diag(ssdragon))

        # Compute p-value
        pval_precision = self.compute_precision_pvalue_normal(ssprecision, df)
        # Compute FDR-corrected p-value using Benjamini-Hochberg
        pval_precision_adjusted = self.benjamini_hochberg(pval_precision)

        if not keep_in_memory:
            print('Saving SMAUG for sample %s' % (str(sample)))
            ssdragon = pd.DataFrame(ssdragon)
            ssprecision = pd.DataFrame(ssprecision)
            pval_precision = pd.DataFrame(pval_precision)
            pval_precision_adjusted = pd.DataFrame(pval_precision_adjusted)
            sfolder = self.output_folder + './smaug/'
            pfolder = self.output_folder + './precision/'
            pvalfolder = self.output_folder + './pval/'
            adjPvalfolder = self.output_folder + './adjPval/'
            if not os.path.exists(sfolder):
                os.makedirs(sfolder)
            if not os.path.exists(pfolder):
                os.makedirs(pfolder)
            if not os.path.exists(pvalfolder):
                os.makedirs(pvalfolder)
            if not os.path.exists(adjPvalfolder):
                os.makedirs(adjPvalfolder)

            output_fn_smaug = sfolder + 'smaug_' + str(sample) + self.output_fmt
            output_fn_precision = pfolder + 'precision_' + str(sample) + self.output_fmt
            output_fn_pval = pvalfolder + 'pval_' + str(sample) + self.output_fmt
            output_fn_adjPval = adjPvalfolder + 'adjPval_' + str(sample) + self.output_fmt
            if self.output_fmt == '.h5':
                ssdragon.to_hdf(output_fn_smaug, key='smaug', index=False)
                ssprecision.to_hdf(output_fn_precision, key='precision', index=False)
                pval_precision.to_hdf(output_fn_pval, key='pval', index=False)
                pval_precision_adjusted.to_hdf(output_fn_adjPval, key='adjPval', index=False)
            elif self.output_fmt == '.csv':
                ssdragon.to_csv(output_fn_smaug, index=False)
                ssprecision.to_csv(output_fn_precision, index=False)
                pval_precision.to_csv(output_fn_pval, index=False)
                pval_precision_adjusted.to_csv(output_fn_adjPval, index=False)
            elif self.output_fmt == '.txt':
                ssdragon.to_csv(output_fn_smaug, index=False, sep='\t')
                ssprecision.to_csv(output_fn_precision, index=False, sep='\t')
                pval_precision.to_csv(output_fn_pval, index=False, sep='\t')
                pval_precision_adjusted.to_csv(output_fn_adjPval, index=False, sep='\t')

            else:
                print('WARNING: output format (%s) not recognised. We are saving in hdf' % str(self.output_fmt))
                ssdragon.to_hdf(output_fn_smaug, key='smaug', index=False)
                ssprecision.to_hdf(output_fn_precision, key='precision', index=False)
                pval_precision.to_hdf(output_fn_pval, key='pval', index=False)
                pval_precision_adjusted.to_hdf(output_fn_adjPval, key='adjPval', index=False)
        else:
            return ssdragon, ssprecision, pval_precision, pval_precision_adjusted

    def compute_precision_pvalue_normal(self, ssprecision, df):
        """
            Compute element-wise p-value

            Parameters:
            ssprecision (2D array-like): sample-specific precision matrix
            df: degrees of freedom of Wishart distribution

            Returns:
            np.ndarray: p-values matrix, same shape as input.
        """

        # Extract the diagonal elements
        diag_elements = np.diag(ssprecision)

        # Create an outer product of the diagonal elements
        outer_diag_product = np.outer(diag_elements, diag_elements)

        v = np.sqrt(outer_diag_product) / df ** 0.75 # under null, off-diagonals of ssprecision = 0

        v = np.divide(ssprecision , v)

        # Two-sided p-value: 2 * (1 - CDF(|z|)) = 2 * SF(|z|)
        p_val = 2 * (1 - norm.cdf(np.abs(v)))

        return p_val

    def benjamini_hochberg(self, pvals):
        """
        Vectorized Benjamini-Hochberg correction for a 2D matrix of p-values.

        Parameters:
        pvals (2D array-like): Input matrix of p-values.

        Returns:
        np.ndarray: Adjusted p-values matrix, same shape as input.
        """
        shape = pvals.shape  # save original shape
        flat_pvals = pvals.flatten()  # flatten to 1D

        # Apply BH correction
        adj_pvals = fdr(flat_pvals, method='bh')

        # Reshape back to 2D
        return adj_pvals.reshape(shape)









