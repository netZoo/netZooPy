from __future__ import print_function
import math
import time
import pandas as pd
from .timer import Timer
import numpy as np
from netZooPy.panda.panda import Panda
from netZooPy.panda import calculations as calc
from netZooPy.ligress import io
import sys
import os
import pandas as pd
import scipy.stats as stats

def compute_bonobo(expression_matrix, expression_mean, sample_idx, online_coexpression = False, computing = 'cpu', cores = 1, delta = None, compute_sparse = False, confidence = 0.05):
    
    """Compute one bonobo matrix. Takes as input an expression matrix, the row-wise mean of the expression, and the
    index of the sample for which the bonobo is computed (index of the columns).
    
    expression_matrix (numpy array): genes as rows, samples as columns
    expression_mean (numpy array): rowise mean of expression. 
    sample_idx (int): index of the sample of interest
    delta (float, optional): delta value for the computation. Defaults to None, which means that delta is tuned
    
    """
    pval = None
    
    mask_include = [True]*expression_matrix.shape[1]
    mask_include[sample_idx] = False
    
    # Compute covariance matrix from the rest of the data, leaving out sample
    covariance_matrix = np.cov(expression_matrix[:, mask_include])
    
    # Compute posterior weight delta from data
    if (delta==None):
        delta = 1/( 3 + 2 * np.sqrt(covariance_matrix.diagonal()).mean()/covariance_matrix.diagonal().var())
    else:
        assert type(delta)==float
    
    # Compute sample-specific covariance matrix
    sscov = delta * np.outer((expression_matrix-expression_mean)[:, sample_idx], (expression_matrix-expression_mean)[:, sample_idx]) + (1-delta) * covariance_matrix

    # Compute sample-specific coexpression matrix from the sample-specific covariance matrix
    
    sscov = np.array(sscov)
    diag = np.sqrt(np.diag(np.diag(sscov)))
    sds = np.linalg.inv(diag)
    bonobo_matrix = sds @ sscov @ sds
    
    
    if compute_sparse:
        threshold  = stats.norm.ppf(1-(confidence/2))
        g = sscov.shape[1]
        d = g + 1/delta
        a1 = (d-g+1)/((d-g)*(d-g-3))
        a2 = (d-g-1)/((d-g)*(d-g-3))
        
        #m1 = a1*(np.multiply(sscov, sscov)) 
        v = np.sqrt(np.diag(sscov)) 
        #m2 = a2 * (np.outer(v, v))
        
        # Pointwise variance: v = m1 + m2
        v = (a1*(np.multiply(sscov, sscov))) + (a2 * (np.outer(v, v)))
        # Pointwise standard deviation
        v = np.sqrt(v)
        
        # zscore
        v = np.divide(sscov,v)
        
        pval = 2*(1-stats.norm.cdf(np.abs(v)))
        
        # bonobo gets sparsified: diagonal bonobo + sparsified off diagonal
        # bonobo_matrix = (np.eye(g)@bonobo_matrix) + np.multiply( 1-np.eye(g), (v*(v<threshold) ) )
        bonobo_matrix = np.eye(g) + np.multiply( 1-np.eye(g), (bonobo_matrix*(np.abs(v)<threshold) ) )

    return(bonobo_matrix, delta, pval)

class Bonobo():
    """
    BONOBO


    Parameters
    ----------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.
            priors_table_file : str
                Path to file containing a table where each samples is linked to its own motif prior file
            ppi_file : str
                Path to file containing the PPI data. or pandas dataframe.
                The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.
            mode_process : str
                The input data processing mode.
                - 'legacy': refers to the processing mode in netZooPy<=0.5
                - (Default)'union': takes the union of all TFs and genes across priors and fills the missing genes in the priors with zeros.
                - 'intersection': intersects the input genes and TFs across priors and removes the missing TFs/genes.
            mode_priors: str
                The prior data processing
            prior_tf_col: str
                name of the tf column in the prior files
            prior_gene_col: str
                name of the gene column in the prior files
            output_folder: str
                folder where to save the results
            delta: float
                posterior weight between 0 and 1 (Default to 0.3)
            tune_delta: boolean
                if true, the posterior weight (delta) for the estimation of the single sample coexpression is estimated from data


    Notes
    ------

    Toy data:The example gene expression data that we have available here contains gene expression profiles
    for different samples in the columns. Of note, this is just a small subset of a larger gene
    expression dataset. We provided these "toy" data so that the user can test the method.


    Sample PANDA results:\b
        - TF    Gene    Motif   Force\n
        - CEBPA AACSL	0.0	-0.951416589143\n
        - CREB1 AACSL	0.0	-0.904241609324\n
        - DDIT3 AACSL	0.0	-0.956471642313\n
        - E2F1  AACSL	1.0	3.685316051\n
        - EGR1  AACSL	0.0	-0.695698519643

    References
    ----------
    .. [1]__ 

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
        print('BONOBO: preparing expression')
        self._prepare_expression()
        self.delta = {}
        self.bonobos = {}
        

    ########################
    ### METHODS ############
    ########################
    
    def _prepare_expression(self):

        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file, samples=self.samples
            )

            self.expression_samples = self.expression_data.columns.tolist()
            
        # Auxiliary dicts
        self.expression_gene2idx = {x: i for i, x in enumerate(self.expression_genes)}

    def run_bonobo(self, output_folder = 'bonobo/', output_fmt = 'hd5', keep_in_memory = False,save_full = False, online_coexpression = False,delta = None, computing = 'cpu', cores = 1, precision = 'single', sample_names = [], sparsify = False, confidence = 0.05):
        
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
        """

        ligress_start = time.time()
        # first let's reorder the expression data
        
        if precision=='single':
            atype = 'float32'
        elif precision=='double':
            atype = 'float64'
        else: 
            sys.exit('Precision %s unknonw' %str(precision))
        
        # let's sort the expression and ppi data
        self.expression_data = self.expression_data.astype(atype)
        
        
        self.sparsify = sparsify
        self.confidence = confidence
        # If output folder is an empty string, keep the matrix in memory and don't save it to disk
        # Otherwise the output folder can be created and the matrix saved
        if output_folder=='':
            keep_in_memory = True
            save_matrix = False
        else:
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            save_matrix = True
        
        
        #correlation_complete = self.expression_data.T.corr()
        # we automatically multiply the correlation with the number of samples
        
        # Center expression data to make mean = 0
        # let's remove this from here, and keep it only inside the ligress computation
        #self.expression_data_centered = (self.expression_data - np.mean(self.expression_data.values,axis = 1, keepdims=True))
        self.expression_mean = np.mean(self.expression_data.values,axis = 1, keepdims=True)
        
        print('BONOBO: We are starting to compute the bonobos...')
        # Now for each sample we compute the lioness network from correlations and 
        # the panda using the motif and ppi tables
        
        if sample_names==[]:
            sample_names=self.expression_samples
        else:
            different = set(sample_names).difference(set(self.expression_samples))
            sample_names = set(sample_names).intersection(set(self.expression_samples))
            if len(different)>0:
                print('WARNING: some of the sample names are not in the expression data')
                print('\tMissing:')
                print('\t'+str(different))
                print('\tUsing:')
                print('\t'+str(sample_names))
        
        
        for s,sample in enumerate(sample_names):
            sample_start = time.time()
            # first run bonobo
            print('BONOBO: bonobo for sample %s' %str(sample))
            self._bonobo_loop(sample, output_fmt = output_fmt, keep_in_memory=keep_in_memory, save_matrix=save_matrix, computing=computing, output_folder=output_folder, delta = delta)

    def _bonobo_loop(self, sample, output_fmt = '.h5', keep_in_memory = False, save_matrix = True, online_coexpression = False, computing= 'cpu', output_folder = './bonobo/' , delta=None):
        """Runs ligress on one sample. For now all samples are saved separately.

        Args:
            correlation_complete (_type_): _description_
            ppi_data (_type_): _description_
            motif_data (_type_): _description_
            sample (_type_): _description_
            keep_coexpression (bool, optional): _description_. Defaults to False.
            save_memory (bool, optional): _description_. Defaults to True.
            online_coexpression (bool, optional): _description_. Defaults to False.
            computing_lioness (str, optional): _description_. Defaults to 'cpu'.
            coexpression_folder (str, optional): _description_. Defaults to './coexpression/'.
            computing_panda (str, optional): _description_. Defaults to 'cpu'.
            alpha (float, optional): _description_. Defaults to 0.1.
        """
        
        touse = list(set(self.expression_samples).difference(set([sample])))
        sample_idx = list(self.expression_samples).index(sample)

        print('BONOBO: computing bonobo for sample %s' %str(sample))
        sample_bonobo, sample_delta, pval = compute_bonobo(self.expression_data.values,self.expression_mean, sample_idx, delta = delta, online_coexpression = online_coexpression, computing = computing, compute_sparse = self.sparsify, confidence = self.confidence)
        
        self.pvals = pval
        
        self.delta[sample] = sample_delta
        
        df_bonobo = pd.DataFrame(data = sample_bonobo, columns = list(self.expression_genes))
        
        if (save_matrix):
            print('Saving BONOBO for sample %s' %(str(sample)))
            output_fn = output_folder + 'bonobo_' + str(sample) + output_fmt
            if output_fmt=='.h5':
                df_bonobo.to_hdf(output_fn,index = False)
            elif output_fmt=='.csv':
                df_bonobo.to_csv(output_fn, index = False)
            elif output_fmt=='.txt':
                df_bonobo.to_csv(output_fn, index = False, sep = '\t')
            else:
                print('WARNING: output format (%s) not recognised. We are saving in hdf' %str(output_fmt))
                output_fn = output_folder + 'bonobo_' + str(sample) + '.h5'
                df_bonobo.to_hdf(output_fn,key ='bonobo',index = False)
        
        if keep_in_memory:
            self.bonobos[sample] = df_bonobo
