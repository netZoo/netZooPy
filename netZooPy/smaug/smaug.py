from __future__ import print_function
import time
from netZooPy.smaug.timer import Timer
import numpy as np
import sys
import os
import pandas as pd
from netZooPy.dragon import estimate_penalty_parameters_dragon
from netZooPy.dragon import get_shrunken_covariance_dragon
from netZooPy.smaug import io


def compute_smaug(expression_matrix, methylation_matrix, expression_mean, methylation_mean, sample_idx, online_partial_coexpression=False, computing='cpu', cores=1,
                   delta=None, compute_sparse=False, confidence=0.05, save_pvals=False):
    """Compute one smaug matrix. Takes as input an expression matrix, a methylation matrix, and the
    index of the sample for which smaug is computed (index of the columns).

    expression_matrix (numpy array): genes as rows, samples as columns
    expression_mean (numpy array): row-wise mean of expression.
    methylation_matrix (numpy array): methylation probes as rows, samples as columns
    methylation_mean (numpy array): row-wise mean of methylation.
    sample_idx (int): index of the sample of interest
    delta (float, optional): delta value for the computation. Defaults to None, which means that delta is tuned

    """

    mask_include = [True] * expression_matrix.shape[1]
    mask_include[sample_idx] = False

    # Compute covariance matrix from the rest of the data, leaving out sample
    lambdas = estimate_penalty_parameters_dragon(expression_matrix[:, mask_include], methylation_matrix[:, mask_include])
    covariance_matrix = get_shrunken_covariance_dragon(expression_matrix[:, mask_include], methylation_matrix[:, mask_include], lambdas)

    # Compute posterior weight delta from data
    if (delta == None):
        delta = 1 / (3 + 2 * np.sqrt(covariance_matrix.diagonal()).mean() / covariance_matrix.diagonal().var())
    else:
        assert type(delta) == float

    # Append expression amd methylation matrix by row
    combined_data = np.concatenate((expression_matrix, methylation_matrix), axis=0)

    # Append mean expression and mean methylation
    combined_mean = np.concatenate((expression_mean, methylation_mean))

    # Compute sample-specific covariance matrix
    sscov = delta * np.outer((combined_data - combined_mean)[:, sample_idx],
                             (combined_data - combined_mean)[:, sample_idx]) + (1 - delta) * covariance_matrix

    # Compute sample-specific DRAGON from sample-specific covariance
    Theta = np.linalg.inv(sscov)
    p = Theta.shape[0]
    A = np.sqrt(np.zeros((p, p)) + np.diag(Theta))
    smaug_matrix = -Theta / A / A.T
    smaug_matrix = smaug_matrix - np.diag(np.diag(smaug_matrix))

    return (smaug_matrix)

class Bonobo():
    """
    BONOBO


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
        self.smaugs = {}
        self.pvals = {}
        self.save_pvals = False

    ########################
    ### METHODS ############
    ########################
    def _prepare_data(self):
        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file, samples=self.samples
            )

        with Timer("Reading methylation data..."):
            # Read expression
            self.methylation_data_data, self.methylation_probes = io.prepare_expression(
                self.methylation_file, samples=self.samples
            )

            self.expression_samples = self.expression_data.columns.tolist()
            self.methylation_samples = self.methylation_data.columns.tolist()

    def run_smaug(self, output_folder='smaug/', output_fmt='hd5', keep_in_memory=False, save_full=False,
                   online_partial_coexpression=False, delta=None, computing='cpu', cores=1, precision='single', sample_names=[],
                   sparsify=False, confidence=0.05, save_pvals=False):
        """BONOBO algorithm

        Args:
            output_folder (str, optional): output folder. If an empty string is passed the matrix is automatically kept
            in memory, overwriting the value of keep_in_memory
            output_fmt (str, optional): format of output bonobo matrix. By default it is an hd5 file, can be a txt or csv.
            keep_in_memory (bool, optional): if True, the partial correlation matrix is kept in memory, otherwise it is
            discarded after saving.
            save_full (bool, optional): whether to save the partial coexpression with the gene names. We recommend using True
            only when the number of genes is not very big to avoid saving huge matrices.
            online_partial_coexpression (bool, optional): if true partial coexpression is computed with a closed form
            cores (int, optional): cores. Defaults to 1.
            delta (float, optional): delta parameter. If default (None) delta is trained, otherwise pass a value.
            precision (str, optional): matrix precision, defaults to single precision.
            sparsify (bool, optional): if True, smaug gets sparsified and relative pvalues are returned
            confidence (float, optional): if sparsify is True, this is the CI for the approximate zscore.
        """

        smaug_start = time.time()

        # first let's reorder the expression data

        if precision == 'single':
            atype = 'float32'
        elif precision == 'double':
            atype = 'float64'
        else:
            sys.exit('Precision %s unknonwn' % str(precision))

        # let's sort the expression and ppi data
        self.expression_data = self.expression_data.astype(atype)
        self.sparsify = sparsify
        self.confidence = confidence
        self.save_pvals = save_pvals
        # If output folder is an empty string, keep the matrix in memory and don't save it to disk
        # Otherwise the output folder can be created and the matrix saved
        if output_folder == '':
            keep_in_memory = True
            save_matrix = False
        else:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            save_matrix = True

        # Compute mean expression and mean methylation
        self.expression_mean = np.mean(self.expression_data.values, axis=1, keepdims=True)
        self.methylation_mean = np.mean(self.methylation_data.values, axis=1, keepdims=True)

        print('SMAUG: We are starting to compute the networks...')
        if sample_names == []:
            sample_names = self.expression_samples
            sample_names = set(sample_names).intersection(set(self.methylation_samples))
        else:
            different = set(sample_names).difference(set(self.expression_samples))
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
            # first run bonobo
            print('SMAUG: network for sample %s' % str(sample))
            self._smaug_loop(sample, output_fmt=output_fmt, keep_in_memory=keep_in_memory, save_matrix=save_matrix,
                              computing=computing, output_folder=output_folder, delta=delta)

    def _bonobo_loop(self, sample, output_fmt='.h5', keep_in_memory=False, save_matrix=True, online_partial_coexpression=False,
                     computing='cpu', output_folder='./smaug/', delta=None):
        """Runs smaug on one sample. All samples are saved separately.

        Args:
            sample (_type_): _description_
            online_partial_coexpression (bool, optional): _description_. Defaults to False.
            computing (str, optional): _description_. Defaults to 'cpu'.
            output_folder (str, optional): _description_. Defaults to './coexpression/'.
        """

        touse = list(set(self.expression_samples).difference(set([sample])))
        sample_idx = list(self.expression_samples).index(sample)

        print('SMAUG: computing network for sample %s' % str(sample))
        sample_smaug, sample_delta, pval = compute_smaug(self.expression_data.values, self.methylation_data.values,
                                                         self.expression_mean, self.methylation_mean,
                                                         sample_idx, delta=delta,
                                                         online_partial_coexpression=online_partial_coexpression,
                                                         computing=computing, compute_sparse=self.sparsify,
                                                         confidence=self.confidence, save_pvals=self.save_pvals)

        self.delta[sample] = sample_delta

        df_smaug = pd.DataFrame(data=sample_smaug, columns=self.expression_data.index.tolist())

        if save_matrix:
            print('Saving SMAUG for sample %s' % (str(sample)))
            output_fn = output_folder + 'smaug_' + str(sample) + output_fmt
            if output_fmt == '.h5':
                df_smaug.to_hdf(output_fn, key='bonobo', index=False)
            elif output_fmt == '.csv':
                df_smaug.to_csv(output_fn, index=False)
            elif output_fmt == '.txt':
                df_smaug.to_csv(output_fn, index=False, sep='\t')
            else:
                print('WARNING: output format (%s) not recognised. We are saving in hdf' % str(output_fmt))
                output_fn = output_folder + 'smaug_' + str(sample) + '.h5'
                df_smaug.to_hdf(output_fn, key='smaug', index=False)

        if (self.sparsify and self.save_pvals):
            df_pval = pd.DataFrame(data=pval, columns=self.expression_data.index.tolist())
            print('Saving pvalues for sample %s' % (str(sample)))
            output_fn = output_folder + 'pvals_' + str(sample) + output_fmt
            if output_fmt == '.h5':
                df_pval.to_hdf(output_fn, key='pvals', index=False)
            elif output_fmt == '.csv':
                df_pval.to_csv(output_fn, index=False)
            elif output_fmt == '.txt':
                df_pval.to_csv(output_fn, index=False, sep='\t')
            else:
                print('WARNING: output format (%s) not recognised. We are saving in hdf' % str(output_fmt))
                output_fn = output_folder + 'pvals_' + str(sample) + '.h5'
                df_pval.to_hdf(output_fn, key='pvals', index=False)

        if keep_in_memory:
            self.smaugs[sample] = df_smaug
            self.pvals[sample] = df_pval






