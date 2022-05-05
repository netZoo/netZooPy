from __future__ import print_function
import math
import time
from tkinter.messagebox import RETRY
import pandas as pd
from scipy.stats import zscore
from .timer import Timer
import numpy as np
from netZooPy.panda.panda import Panda
from netZooPy.panda import calculations as calc
from netZooPy.tigress import io
import sys
import os


class Tigress(Panda):
    """
    Using PANDA to infer gene regulatory network from coexpression networks.
        1. Reading in input data (coexpression, motif prior, TF PPI data)
        3. Normalizing networks
        4. Running PANDA algorithm
        5. Writing out PANDA network


    Parameters
    ----------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.
            motif_file : str
                Path to file containing the transcription factor DNA binding motif data in the form of
                TF-gene-weight(0/1) or pandas dataframe.
                If set to none, the gene coexpression matrix is returned as a result network.
            ppi_file : str
                Path to file containing the PPI data. or pandas dataframe.
                The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.
            computing : str
                'cpu' uses Central Processing Unit (CPU) to run PANDA.
                'gpu' use the Graphical Processing Unit (GPU) to run PANDA.
            precision : str
                - 'double' computes the regulatory network in double precision (15 decimal digits).
                - 'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate.
            save_memory : bool
                - True : removes temporary results from memory. The result network is weighted adjacency matrix of size (nTFs, nGenes).
                - False: keeps the temporary files in memory. The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.
            save_tmp : bool
                Save temporary variables.
            remove_missing : bool
                Removes the gens and TFs that are not present in one of the priors. Works only if modeProcess='legacy'.
            keep_expression_matrix : bool
                Keeps the input expression matrix in the result Panda object.
            modeProcess : str
                The input data processing mode.
                - 'legacy': refers to the processing mode in netZooPy<=0.5
                - (Default)'union': takes the union of all TFs and genes across priors and fills the missing genes in the priors with zeros.
                - 'intersection': intersects the input genes and TFs across priors and removes the missing TFs/genes.
            alpha : str
                Learning rate (default: 0.1)
            start : int
                First sample of the expression dataset. This replicates the behavior of Lioness (default : 1)
            end : int
            Last sample of the expression dataset. This replicates the behavior of Lioness (default : None )

    Examples
    --------
        # TODO: fix here
        >>> #Import the classes in the pypanda library:
        >>> from netZooPy.panda.panda import Panda
        >>> #Run the Panda algorithm, leave out motif and PPI data to use Pearson correlation network:
        >>> panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False)
        >>> #Save the results:
        >>> panda_obj.save_panda_results('Toy_Panda.pairs.txt')
        >>> #Return a network plot:
        >>> panda_obj.top_network_plot(top=70, file='top_genes.png')
        >>> #Calculate in- and outdegrees for further analysis:
        >>> indegree = panda_obj.return_panda_indegree()
        >>> outdegree = panda_obj.return_panda_outdegree()


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
    .. [1]__ Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions."
        PloS one 8.5 (2013): e64832.

    Authors: Cho-Yi Chen, David Vi, Alessandro Marin, Marouen Ben Guebila, Daniel Morgan

    """

    def __init__(
        self,
        expression_file,
        priors_table_file,
        ppi_file,
        mode_process="union",
        mode_priors="union",
        prior_tf_col=0,
        prior_gene_col=1,
        output_folder='./tigress/' 
    ):
        """Intialize instance of Panda class and load data."""

        self.expression_file = expression_file
        self.priors_table_file = priors_table_file
        self.ppi_file = ppi_file
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
            # prepare universe of names in the priors. We won't be reading all of them 
            # first, because we might want to use too many motif priors
            (
                self.priors_tfs,
                self.priors_genes,
            ) = io.read_motif_universe(
                self.sample2prior_dict, mode=self.mode_priors
            )

            

        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file, self.samples
            )

        with Timer("Loading PPI data ..."):
            # read ppi
            self.ppi_data, self.ppi_tfs = io.read_ppi(self.ppi_file)

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
        self.expression_data = self.expression_data.loc[self.universe_genes,:]
        self.ppi_data = self.ppi_data.loc[self.universe_tfs,self.universe_tfs]
    
    def run_tigress(self, keep_coexpression = False,save_memory = False, online_coexpression = False, coexpression_folder = 'coexpression/', computing_lioness = 'cpu', computing_panda = 'cpu', cores = 1, alpha = 0.1 , precision = 'single', th_motifs = 3):
        
        """Tigress algorithm

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

        tigress_start = time.time()
        # first let's reorder the expression data
        
        if precision=='single':
            atype = 'float32'
        elif precision=='double':
            atype = 'float64'
        else: 
            sys.exit('Precision %s unknonw' %str(precision))
        
        # let's sort the expression and ppi data

        self.expression_data = self.expression_data.loc[self.universe_genes,:].astype(atype)
        correlation_complete = self.expression_data.T.corr()
        

        if th_motifs>len(self.prior2sample_dict.keys()):
            for p,ss in self.prior2sample_dict.items():
                # read the motif data and sort it
                motif_data = self._get_motif(p)
                for s,sample in enumerate(ss):
                    sample_start = time.time()
                    # first run lioness on coexpression
                    self._tigress_loop(correlation_complete, motif_data, sample, keep_coexpression=keep_coexpression, save_memory=save_memory, computing_lioness=computing_lioness, computing_panda=computing_panda, alpha = alpha, coexpression_folder=coexpression_folder)

        else:
            # Now for each sample we compute the lioness network from correlations and 
            # the panda using the motif and ppi tables
            for s,sample in enumerate(self.samples):
                sample_start = time.time()
                # first run lioness on coexpression
                motif_data = self._get_motif(self.sample2prior_dict[sample])
                self._tigress_loop(correlation_complete, motif_data, sample, keep_coexpression=keep_coexpression, save_memory=save_memory, computing_lioness=computing_lioness, computing_panda=computing_panda, alpha = alpha, coexpression_folder=coexpression_folder)


    def _tigress_loop(self, correlation_complete, motif_data, sample, keep_coexpression = False, save_memory = True, online_coexpression = False, computing_lioness = 'cpu', coexpression_folder = './coexpression/' , computing_panda = 'cpu', alpha = 0.1):
        
        if keep_coexpression:
            if not os.path.exists(self.output_folder+coexpression_folder):
                os.makedirs(self.output_folder+coexpression_folder)
        
        if not os.path.exists(self.output_folder+'single_panda/'):
            os.makedirs(self.output_folder+'single_panda/')
        sample_lioness = self._run_lioness_coexpression(correlation_complete, sample, keep_coexpression = keep_coexpression, save_memory = save_memory, online_coexpression = online_coexpression, computing = computing_lioness, coexpression_folder = coexpression_folder)

        final_panda= self._run_panda_coexpression(sample_lioness,motif_data, sample, computing = computing_panda, alpha = alpha, save_single=True)

    def _save_single_panda_net(self, net, prior, sample, prefix, pivot = False):

        tab = pd.DataFrame(net, columns = self.universe_genes )
        tab['tf'] = self.universe_tfs

        if pivot:
            tab.set_index('tf').to_csv(prefix+sample+'.csv')
        else:
            tab = pd.melt(tab, id_vars='tf', value_vars=tab.columns,var_name='gene', value_name='force')
            tab['motif'] = prior.flatten()
            tab.to_csv(prefix+sample+'.txt', sep = '\t', index = False, columns = ['tf', 'gene','motif','force'])

    def _get_motif(self, motif_fn):
        motif_data = io.read_motif(motif_fn, tf_names = list(self.universe_tfs), gene_names = list(self.universe_genes), pivot = True)
        return(motif_data)

    def _run_lioness_coexpression(self, coexpression, sample, keep_coexpression = False,save_memory = True, online_coexpression = False, computing = 'cpu', cores = 1, coexpression_folder = 'coexpression/'):
        
        touse = set(self.samples).difference(set([sample]))
        names = self.expression_data.index.tolist()
        correlation_matrix = self.expression_data.loc[:, touse].T.corr().values

        
        
        # For consistency with R, we are using the N panda_all - (N-1) panda_all_but_q
        lioness_network = (self.n_samples * coexpression) - (
                (self.n_samples - 1) * correlation_matrix
        )

        if (keep_coexpression):
            cfolder = self.output_folder+coexpression_folder
            if not os.path.exists(cfolder):
                os.makedirs(cfolder)
            if (save_memory):
                np.savetxt(coexpression_folder+'coexpression_'+sample+'.txt', lioness_network)
            else:
                pd.DataFrame(lioness_network, columns=self.universe_genes, index = self.universe_genes).to_csv(cfolder+'coexpression_'+sample+'.txt', sep = ' ')
        
        return(lioness_network)

    def _run_panda_coexpression(self, net, motif, sample, computing = 'cpu', alpha = 0.1, save_single = False):
        
        panda_loop_time = time.time()
        
        #panda works with all normalised networks
        final = calc.compute_panda(
            self._normalize_network(net.values),
            self._normalize_network(self.ppi_data.astype(float).values),
            self._normalize_network(motif.astype(float).values),
            computing=computing,
            alpha=alpha,
        )
        print("Running panda took: %.2f seconds!" % (time.time() - panda_loop_time))

        if save_single:

            self._save_single_panda_net(final, motif.values, sample, prefix = self.output_folder+'single_panda/', pivot = False)
        return(final)

