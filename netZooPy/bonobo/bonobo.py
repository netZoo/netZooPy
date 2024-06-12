from __future__ import print_function
import math
import time
import pandas as pd
from .timer import Timer
import numpy as np
from netZooPy.panda.panda import Panda
from netZooPy.panda import calculations as calc
from netZooPy.bonobo import io
import sys
import os
import pandas as pd
import scipy.stats as stats


def compute_bonobo(
    expression_matrix,
    expression_mean,
    sample_idx,
    online_coexpression=False,
    computing="cpu",
    cores=1,
    delta=None,
    compute_sparse=False,
    confidence=0.05,
    save_pvals=False,
):

    """Compute one bonobo matrix. Takes as input an expression matrix, the row-wise mean of the expression, and the
    index of the sample for which the bonobo is computed (index of the columns).

    Args:
        expression_matrix (numpy array): genes as rows, samples as columns  
        expression_mean (numpy array): rowise mean of expression.  
        sample_idx (int): index of the sample of interest
        delta (float, optional): delta value for the computation. Defaults to None, which means that delta is tuned  
        compute_sparse (bool, optional): if True, the bonobo gets sparsified and relative pvalues are returned  
        confidence (float, optional): if sparsify is True, this is the CI for the approximate zscore.  
        save_pvals (bool, optional): if True, the pvalues are saved and returned  
        online_coexpression (bool, optional): FOR FUTURE use. if True, the coexpression is computed with a closed form  
        computing (str, optional): FOR FUTURE use. Defaults to 'cpu'.  
        cores (int, optional): FOR FUTURE use. Defaults to 1.    
          
    Return:
        bonobo_matrix: bonobo coexpression matrix  
        delta: delta value, useful when tuned  
        pval: pvalues, if save_pvals is True  
    """
    
    pval = None

    mask_include = [True] * expression_matrix.shape[1]
    mask_include[sample_idx] = False

    # Compute covariance matrix from the rest of the data, leaving out sample
    covariance_matrix = np.cov(expression_matrix[:, mask_include])

    # Compute posterior weight delta from data
    if delta == None:
        delta = 1 / (
            3
            + 2
            * np.sqrt(covariance_matrix.diagonal()).mean()
            / covariance_matrix.diagonal().var()
        )
    else:
        assert type(delta) == float

    # Compute sample-specific covariance matrix
    sscov = (
        delta
        * np.outer(
            (expression_matrix - expression_mean)[:, sample_idx],
            (expression_matrix - expression_mean)[:, sample_idx],
        )
        + (1 - delta) * covariance_matrix
    )

    # Compute sample-specific coexpression matrix from the sample-specific covariance matrix

    sscov = np.array(sscov)
    diag = np.sqrt(np.diag(np.diag(sscov)))

    # Replace 0 diagonals by 1, so that the diagonal matrix can be inverted
    diag = np.array(diag)
    indices = np.where(np.diag(diag) == 0)[0]
    for i in indices:
        diag[i, i] = 1

    sds = np.linalg.inv(diag)
    bonobo_matrix = sds @ sscov @ sds

    if compute_sparse:
        threshold = stats.norm.ppf(1 - (confidence / 2))
        g = sscov.shape[1]
        d = g + 1 / delta
        a1 = (d - g + 1) / ((d - g) * (d - g - 3))
        a2 = (d - g - 1) / ((d - g) * (d - g - 3))

        # m1 = a1*(np.multiply(sscov, sscov))
        v = np.diag(sscov)
        # m2 = a2 * (np.outer(v, v))

        # Pointwise variance: v = m1 + m2
        v = (a1 * (np.multiply(sscov, sscov))) + (a2 * (np.outer(v, v)))
        # Pointwise standard deviation
        v = np.sqrt(v)

        # zscore
        v = np.divide(sscov, v)

        pval = 2 * (1 - stats.norm.cdf(np.abs(v)))

        # bonobo gets sparsified: diagonal bonobo + sparsified off diagonal
        # bonobo_matrix = (np.eye(g)@bonobo_matrix) + np.multiply( 1-np.eye(g), (v*(v<threshold) ) )
        if save_pvals:
            print("keep Bonobo whole, and use pvals to threshold and sparsify")
        else:
            bonobo_matrix = np.eye(g) + np.multiply(
                1 - np.eye(g), (bonobo_matrix * (np.abs(v) > threshold))
            )

    return (bonobo_matrix, delta, pval)


class Bonobo:
    """
    BONOBO


    Parameters
    ----------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.


    Notes
    ------

        BONOBO can be run as follows:    
            bonobo_obj_sparse = Bonobo(expression_file)
            bonobo_obj_sparse.run_bonobo(keep_in_memory=True, output_fmt='.hdf', sparsify=True, output_folder='../data/processed/bonobo_sparse_pvals/', save_pvals=False)



    References
    ----------
    .. [1]__ Bayesian Optimized sample-specific Networks Obtained By Omics data (BONOBO), Saha and Fanfani et al. 2021. doi: https://doi.org/10.1101/2023.11.16.567119

    Authors: Viola Fanfani, Enakshi Saha
    """

    def __init__(
        self,
        expression_file,
    ):
        """Intialize instance of Panda class and load data."""

        self.expression_file = expression_file

        # data read
        self.samples = None
        self.n_samples = None
        self.expression_data = None
        self.expression_genes = None
        self.expression_gene2idx = None
        self.expression_samples = None
        # prepare all the data
        print("BONOBO: preparing expression")
        self._prepare_expression()
        self.delta = {}
        self.bonobos = {}
        self.pvals = {}
        self.save_pvals = False

    ########################
    ### METHODS ############
    ########################

    def _prepare_expression(self):
        
        """Prepare expression data. Uses the io module to read the expression data."""

        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file, samples=self.samples
            )

            self.expression_samples = self.expression_data.columns.tolist()

        # Auxiliary dicts
        self.expression_gene2idx = {x: i for i, x in enumerate(self.expression_genes)}

    def run_bonobo(
        self,
        output_folder="bonobo/",
        output_fmt=".h5",
        keep_in_memory=False,
        save_full=False,
        online_coexpression=False,
        delta=None,
        computing="cpu",
        cores=1,
        precision="single",
        sample_names=[],
        sparsify=False,
        confidence=0.05,
        save_pvals=False,
    ):

        """BONOBO algorithm

        Args:
            output_folder (str, optional): output folder. If an empty string is passed the matrix is automatically kept
            in memory, overwriting the value of keep_in_memory
            output_fmt (str, optional): format of output bonobo matrix. By default it is an hd5 file, can be a txt or csv.
            keep_in_memory (bool, optional): if True, the bonobo coexpression matrix is kept in memory, otherwise it is
            discarded after saving.
            save_full (bool, optional): whether to save the coexpression with the gene names. We recommend using True
            only when the number of genes is not very big to avoid saving huge matrices.
            online_coexpression (bool, optional): if true coexpression is computed with a closed form
            cores (int, optional): cores. Defaults to 1.
            delta (float, optional): delta parameter. If default (None) delta is trained, otherwise pass a value.
            Recommended is 0.3.
            precision (str, optional): matrix precision, defaults to single precision.
            sparsify (bool, optiona): if True, bonobo gets sparsified and relative pvalues are returned
            confidence (float, optional): if sparsify is True, this is the CI for the approximate zscore.
            save_pvals (bool, optional): if True, the pvalues are saved and returned
        """

        bonobo_start = time.time()

        # first let's reorder the expression data

        if precision == "single":
            atype = "float32"
        elif precision == "double":
            atype = "float64"
        else:
            sys.exit("ERROR: Precision %s unknonw" % str(precision))

        # let's sort the expression and ppi data
        self.expression_data = self.expression_data.astype(atype)

        self.sparsify = sparsify
        self.confidence = confidence
        self.save_pvals = save_pvals
        # If output folder is an empty string, keep the matrix in memory and don't save it to disk
        # Otherwise the output folder can be created and the matrix saved
        if output_folder == "":
            keep_in_memory = True
            save_matrix = False
        else:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            save_matrix = True

        # correlation_complete = self.expression_data.T.corr()
        # we automatically multiply the correlation with the number of samples

        # Center expression data to make mean = 0
        # let's remove this from here, and keep it only inside the bonobo computation
        # self.expression_data_centered = (self.expression_data - np.mean(self.expression_data.values,axis = 1, keepdims=True))
        self.expression_mean = np.mean(
            self.expression_data.values, axis=1, keepdims=True
        )

        print("BONOBO: We are starting to compute the bonobos...")
        # Now for each sample we compute the lioness network from correlations and
        # the panda using the motif and ppi tables

        if sample_names == []:
            sample_names = self.expression_samples
        else:
            different = set(sample_names).difference(set(self.expression_samples))
            sample_names = set(sample_names).intersection(set(self.expression_samples))
            if len(different) > 0:
                print(
                    "WARNING: some of the sample names are not in the expression data"
                )
                print("\tMissing:")
                print("\t" + str(different))
                print("\tUsing:")
                print("\t" + str(sample_names))

        for s, sample in enumerate(sample_names):
            sample_start = time.time()
            # first run bonobo
            print("BONOBO: bonobo for sample %s" % str(sample))
            self._bonobo_loop(
                sample,
                output_fmt=output_fmt,
                keep_in_memory=keep_in_memory,
                save_matrix=save_matrix,
                computing=computing,
                output_folder=output_folder,
                delta=delta,
            )

    def _bonobo_loop(
        self,
        sample,
        output_fmt=".h5",
        keep_in_memory=False,
        save_matrix=True,
        online_coexpression=False,
        computing="cpu",
        output_folder="./bonobo/",
        delta=None,
    ):
        """Runs BONOBO on one sample. For now all samples are saved separately.

        Args:
            sample (str): sample of interest
            output_fmt (str, optional): format of output bonobo matrix. By default it is an hd5 file, can be a txt or csv. Defaults to '.h5'.
            keep_in_memory (bool, optional): If true, all BONOBOs are kept in memory and accessed from the bonobo object. Defaults to False.
            save_matrix (bool, optional): If true, the BONOBO is saved to disk. Defaults to True.
            output_folder (str, optional): results folder. Defaults to './bonobo/'.
            delta (float, optional): Delta value, if None the value is computed with the optimization strategy. Defaults
            to None.
            computing (str, optional): For the future. Defaults to 'cpu'.
            online_coexpression (bool, optional): For the future. Defaults to False.
            
        """

        touse = list(set(self.expression_samples).difference(set([sample])))
        sample_idx = list(self.expression_samples).index(sample)

        print("BONOBO: computing bonobo for sample %s" % str(sample))
        sample_bonobo, sample_delta, pval = compute_bonobo(
            self.expression_data.values,
            self.expression_mean,
            sample_idx,
            delta=delta,
            online_coexpression=online_coexpression,
            computing=computing,
            compute_sparse=self.sparsify,
            confidence=self.confidence,
            save_pvals=self.save_pvals,
        )

        self.delta[sample] = sample_delta

        df_bonobo = pd.DataFrame(
            data=sample_bonobo, columns=self.expression_data.index.tolist()
        )

        if save_matrix:
            print("Saving BONOBO for sample %s" % (str(sample)))
            output_fn = output_folder + "bonobo_" + str(sample) + output_fmt
            if output_fmt == ".h5":
                df_bonobo.to_hdf(output_fn, key="bonobo", index=False)
            elif output_fmt == ".csv":
                df_bonobo.to_csv(output_fn, index=False)
            elif output_fmt == ".txt":
                df_bonobo.to_csv(output_fn, index=False, sep="\t")
            else:
                print(
                    "WARNING: output format (%s) not recognised. We are saving in hdf"
                    % str(output_fmt)
                )
                output_fn = output_folder + "bonobo_" + str(sample) + ".h5"
                df_bonobo.to_hdf(output_fn, key="bonobo", index=False)

        if self.sparsify and self.save_pvals:
            df_pvals = pd.DataFrame(
                data=pval, columns=self.expression_data.index.tolist()
            )
            print("Saving pvalues for sample %s" % (str(sample)))
            output_fn = output_folder + "pvals_" + str(sample) + output_fmt
            if output_fmt == ".h5":
                df_pvals.to_hdf(output_fn, key="pvals", index=False)
            elif output_fmt == ".csv":
                df_pvals.to_csv(output_fn, index=False)
            elif output_fmt == ".txt":
                df_pvals.to_csv(output_fn, index=False, sep="\t")
            else:
                print(
                    "WARNING: output format (%s) not recognised. We are saving in hdf"
                    % str(output_fmt)
                )
                output_fn = output_folder + "pvals_" + str(sample) + ".h5"
                df_pvals.to_hdf(output_fn, key="pvals", index=False)

        if keep_in_memory:
            self.bonobos[sample] = df_bonobo
            self.pvals[sample] = pval
