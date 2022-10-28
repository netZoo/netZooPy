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


class Ligress(Panda):
    """
    Lioness-Inferred Gene REgulatory Sample-Specific networks 
        1. Reading in input data (expression, motif prior table, TF PPI data)
        2. Preparing motif prior universe
        3. Estimating sample-specific coexpression with lioness
        4. Running PANDA with a different prior for each samples
    
    Warning: if you are familiar with the other netzoopy functions, this one is slightly different. 
    We have separated the reading and preprocessing steps from those for computation in a more 
    OOP-friendly fashion.


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
        priors_table_file,
        ppi_table_file=None,
        ppi_file = None,
        mode_process="union",
        mode_priors="union",
        prior_tf_col=0,
        prior_gene_col=1,
        output_folder='./ligress/'
    ):
        """Intialize instance of Panda class and load data."""

        self.expression_file = expression_file
        self.priors_table_file = priors_table_file
        self.ppi_table_file = ppi_table_file
        self.ppi_file = ppi_file
        if self.ppi_file:
            self.ppi_mode = 'motif'
        else:
            self.ppi_mode = 'sample'
        self.mode_process = mode_process
        self.mode_priors = mode_priors
        self.prior_tf_col = prior_tf_col
        self.prior_gene_col=prior_gene_col
        self.output_folder = output_folder


        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        # data read
        self.samples = None
        self.n_samples = None
        self.prior_dict = None
        self.expression_data = None
        self.expression_genes = None

        # we need to keep track of all the names in the expression and motif data
        self.priors_tfs = None
        self.priors_genes = None

        # dictionaries mapping sample:prior_file and prior_file:[samples]
        self.sample2prior_dict = None
        self.prior2sample_dict = None

        # SORTED LIST OF TFS AND GENES
        self.universe_tfs = None
        self.universe_genes = None
        self.gene2idx = None
        self.tf2idx = None


        # prepare all the data
        self._prepare_data()

    ########################
    ### METHODS ############
    ########################
    def _prepare_data(self):

        # Read the sample-prior table. We need to know what samples we are using
        with Timer("Reading sample-prior configuration..."):
            self.samples, self.sample2prior_dict, self.prior2sample_dict = io.read_priors_table(self.priors_table_file)
            self.n_samples = len(self.samples)

            if self.ppi_mode=='sample':
                self.samples_ppi, self.sample2ppi_dict, self.ppi2sample_dict = io.read_priors_table(self.ppi_table_file, sample_col = 'sample', prior_col = 'prior')
                #TODO: add check that samples ppi == samples motif

            # prepare universe of names in the priors. We won't be reading all of them 
            # first, because we might want to use too many motif priors

            # from motif data
            (
                self.priors_tfs,
                self.priors_genes,
            ) = io.read_motif_universe(
                self.sample2prior_dict, mode=self.mode_priors
            )
            
            #from ppi data
            # if ppi table file is specified
            with Timer("Loading PPI data ..."):
                if self.ppi_mode=='sample':
                    (
                        self.ppi_tfs,
                    ) = io.read_ppi_universe(
                        self.sample2ppi_dict, mode=self.mode_priors
                    )
                    self.ppi_data = None
                else:
                    # read ppi
                    self.ppi_data, self.ppi_tfs = io.read_ppi(self.ppi_file)


        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file, samples=self.samples
            )



        # depending on the strategy for 
        if self.mode_process=='intersection':
            self.universe_genes = sorted(list(set(self.expression_genes).intersection(set(self.priors_genes))))
            self.universe_tfs = sorted(list(set(self.ppi_tfs).intersection(set(self.priors_tfs))))
        else:
            sys.exit('Only intersection is an available modeProcess for the moment')

        # Auxiliary dicts
        self.gene2idx = {x: i for i, x in enumerate(self.universe_genes)}
        self.tf2idx = {x: i for i, x in enumerate(self.universe_tfs)}

        # sort the gene expression and ppi data
        self.expression_data = self.expression_data.loc[self.universe_genes,self.samples]

        if self.ppi_mode=='motif':
            self.ppi_data = self.ppi_data.loc[self.universe_tfs,self.universe_tfs]
    
    def run_ligress(self, keep_coexpression = False,save_memory = False, online_coexpression = False, coexpression_folder = 'coexpression/', computing_lioness = 'cpu', computing_panda = 'cpu', cores = 1, alpha = 0.1 , precision = 'single', th_motifs = 3, tune_delta=False, delta=0.1):
        
        """Ligress algorithm

        Args:
            keep_coexpression (bool, optional): whether to save each coexpression network
            save_memory (bool, optional): whether to save the coexpression with the gene names
            online_coexpression (bool, optional): if true coexpression is computed with a closed form
            coexpression_folder (str, optional): used if keep_coexpression is passed
            computing_lioness (str, optional): computing for coexpression lioness. Defaults to 'cpu'.
            computing_panda (str, optional): computing for single sample panda. Defaults to 'cpu'.
            cores (int, optional): cores. Defaults to 1.
            alpha (float, optional): _description_. Defaults to 0.1.
            th_motifs (int, optional): if the number of motif files is lower than the threshold, each will be loaded
            only once.
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

        self.expression_data = self.expression_data.loc[self.universe_genes,:].astype(atype)
        #correlation_complete = self.expression_data.T.corr()
        # we automatically multiply the correlation with the number of samples
        #correlation_complete = correlation_complete * self.get_n_matrix(self.expression_data)
        self.n_matrix = self.get_n_matrix(self.expression_data)
        
        # consider removing this
        # scale expression data to make mean = 0 and sd = 1
        # self.expression_data_scaled = (self.expression_data - self.expression_data.mean(axis = 1))/self.expression_data.std(ddof=1, axis = 1)
    
        # Center expression data to make mean = 0
        # let's remove this from here, and keep it only inside the ligress computation
        #self.expression_data_centered = (self.expression_data - np.mean(self.expression_data.values,axis = 1, keepdims=True))
        self.expression_mean = np.nanmean(self.expression_data.values,axis = 1, keepdims=True)
        
        if th_motifs>len(self.prior2sample_dict.keys()):
            for p,ss in self.prior2sample_dict.items():
                # read the motif data and sort it
                motif_data, tftoadd, genetoadd = self._get_motif(p)
                for s,sample in enumerate(ss):
                    sample_start = time.time()
                    ppi_data = self._get_ppi(sample, missing_tf = tftoadd)
                    # first run lioness on coexpression
                    self._ligress_loop(ppi_data, motif_data, sample, keep_coexpression=keep_coexpression, save_memory=save_memory, computing_lioness=computing_lioness, computing_panda=computing_panda, alpha = alpha, coexpression_folder=coexpression_folder, delta = delta, tune_delta=tune_delta)

        else:
            # Now for each sample we compute the lioness network from correlations and 
            # the panda using the motif and ppi tables
            for s,sample in enumerate(self.samples):
                sample_start = time.time()
                # first run lioness on coexpression
                motif_data, tftoadd, genetoadd = self._get_motif(self.sample2prior_dict[sample])
                ppi_data = self._get_ppi(sample, missing_tf = tftoadd)
                self._ligress_loop(ppi_data, motif_data, sample, keep_coexpression=keep_coexpression, save_memory=save_memory, computing_lioness=computing_lioness, computing_panda=computing_panda, alpha = alpha, coexpression_folder=coexpression_folder, delta = delta, tune_delta=tune_delta)


    def _ligress_loop(self, ppi_data, motif_data, sample, keep_coexpression = False, save_memory = True, online_coexpression = False, computing_lioness = 'cpu', coexpression_folder = './coexpression/' , computing_panda = 'cpu', alpha = 0.1, delta=0.3,tune_delta=False):
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
        if keep_coexpression:
            if not os.path.exists(self.output_folder+coexpression_folder):
                os.makedirs(self.output_folder+coexpression_folder)
        
        if not os.path.exists(self.output_folder+'single_panda/'):
            os.makedirs(self.output_folder+'single_panda/')
        sample_lioness = self._run_lioness_coexpression(sample, keep_coexpression = keep_coexpression, save_memory = save_memory, online_coexpression = online_coexpression, computing = computing_lioness, coexpression_folder = coexpression_folder, delta = delta, tune_delta = tune_delta)

        final_panda= self._run_panda_coexpression(sample_lioness,ppi_data, motif_data, sample, computing = computing_panda, alpha = alpha, save_single=True)
        #return(final_panda)

    def _save_single_panda_net(self, net, prior, sample, prefix, pivot = False):

        tab = pd.DataFrame(net, columns = self.universe_genes )
        tab['tf'] = self.universe_tfs

        if pivot:
            tab.set_index('tf').to_csv(prefix+sample+'.csv')
        else:
            tab = pd.melt(tab, id_vars='tf', value_vars=tab.columns,var_name='gene', value_name='force')
            tab['motif'] = prior.flatten(order = 'F')
            tab.to_csv(prefix+sample+'.txt', sep = '\t', index = False, columns = ['tf', 'gene','motif','force'])

    def _get_motif(self, motif_fn):
        motif_data,tftoadd, genetoadd = io.read_motif(motif_fn, tf_names = list(self.universe_tfs), gene_names = list(self.universe_genes), pivot = True)
        return(motif_data, tftoadd,genetoadd)

    def _get_ppi(self, sample, missing_tf = None):
        if (self.ppi_mode == 'sample'):
            data = io.read_ppi(self.sample2ppi_dict[sample], self.universe_tfs)
        else:
            data = self.ppi_data
            # if there are missing tf, the ppi is all null and 
            if missing_tf:
                data.loc[missing_tf,:]=0
                data.loc[:,missing_tf]=0
                #data.loc[missing_tf,missing_tf] = np.eye(len(missing_tf))

        return(data)

    def _run_lioness_coexpression(self, sample, keep_coexpression = False,save_memory = True, online_coexpression = False, computing = 'cpu', cores = 1, coexpression_folder = 'coexpression/', delta = 0.3, tune_delta = False):
        
        touse = list(set(self.samples).difference(set([sample])))
        names = self.expression_data.index.tolist()
                 
        #correlation_matrix = self.expression_data.loc[:, touse].T.corr().values
        # Compute covariance matrix from the rest of the data, leaving out sample
        covariance_matrix = self.expression_data.loc[:, touse].T.cov().values
        
        # Compute posterior weight delta from data
        if (tune_delta):
            delta = 1/( 3 + 2 * np.nanmean(np.sqrt(covariance_matrix.diagonal()))/np.nanvar(covariance_matrix.diagonal()))
            print(delta)
        # #TODO: remove this old version
        # For consistency with R, we are using the N panda_all - (N-1) panda_all_but_q
        # coexpression has been already multiplied by N all
        # we no longer need coexpression
        #lioness_network = coexpression - (
        #        (self.get_n_matrix(self.expression_data.loc[:, touse])) * correlation_matrix
        #)
        
        # Compute sample-specific covariance matrix
        sscov = delta * np.outer((self.expression_data-self.expression_mean).loc[:, sample], (self.expression_data-self.expression_mean).loc[:, sample]) + (1-delta) * covariance_matrix

        # Compute sample-specific coexpression matrix from the sample-specific covariance matrix
        
        # compute a diagonal matrix from the inverse square root of sscov
        sscov = np.array(sscov)
        print(sscov)
        sscov = np.where(~np.isnan(sscov),sscov,0)
        sscov_diag = np.diag(sscov)
        sds = np.diag(1/np.sqrt(np.where(sscov_diag!=0, sscov_diag, 1)))
        # SS coexpression
        coexp = np.matmul(sds, np.matmul(sscov,sds))
        ## matrix of number of samples (allows to put zeros where NaN)
        nmatrix = self.n_matrix - self.get_n_matrix(self.expression_data.loc[:, touse])
        print(nmatrix)
        coexp = np.multiply(nmatrix, coexp)
        print(coexp)
        #

        if (keep_coexpression):
            cfolder = self.output_folder+coexpression_folder
            if not os.path.exists(cfolder):
                os.makedirs(cfolder)
            path = cfolder+'coexpression_'+sample
            path_genename = cfolder+'genenames_'+sample
            if (save_memory):
                #if self.save_fmt == "txt":
                #np.savetxt(path+'.txt', coexp)
                #elif self.save_fmt == "npy":
                np.save(path+'.npy', coexp)
                # write the gene names
                with open(path_genename+'.txt', 'w') as fp:
                    for item in names:
                        # write each item on a new line
                        fp.write("%s\n" % item)
                #elif self.save_fmt == "mat":
                #    from scipy.io import savemat
                #    savemat(path, {"SSCoexp": coexp})
            else:
                pd.DataFrame(data = coexp, columns=names, index = names).to_csv(cfolder+'coexpression_'+sample+'.txt', sep = ' ')
        
        return(pd.DataFrame(data = coexp, index = names, columns=names))



    def _run_panda_coexpression(self, net, ppi, motif, sample, computing = 'cpu', alpha = 0.1, save_single = False):
        
        panda_loop_time = time.time()
        
        #panda works with all normalised networks
        if (len(ppi.index)!=np.sum(ppi.index==motif.index)):
            sys.exit('PPI and motif tfs are not matching. DEBUG!')
        if (len(net.index)!=np.sum(motif.columns==net.index)):
            sys.exit('coexpression and motif genes are not matching. DEBUG!')
        final = calc.compute_panda(
            self._normalize_network(net.values),
            self._normalize_network(ppi.values),
            self._normalize_network(motif.astype(float).values),
            computing=computing,
            alpha=alpha,
        )
        print("Running panda took: %.2f seconds!" % (time.time() - panda_loop_time))

        if save_single:
            self._save_single_panda_net(final, motif.values, sample, prefix = self.output_folder+'single_panda/', pivot = False)
        return(final)

    def get_n_matrix(self,df):
        # This should be outside of the class
        """Get number of samples for each correlation value

        Args:
            df (pd.DataFrame): expression with nan values
        """
        
        N = len(df.columns)
        nn = N-df.isna().sum(axis = 1).values[:,np.newaxis]
        nr = np.repeat(nn,len(nn), axis = 1)
        return(np.minimum(nr,nr.T))
