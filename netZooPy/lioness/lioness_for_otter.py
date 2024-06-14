from __future__ import print_function
import math
import time
import pandas as pd
from .timer import Timer
import numpy as np
from netZooPy.otter.otter import otter
from netZooPy.otter.otter import otter_gpu
from netZooPy.lioness import io
import sys
import os


class LionessOtter():
    """
    Lioness-Inferred Single-sample Otter GRNs. 
        1. Reading in input data (expression, motif prior table, TF PPI data)
        2. Running Otter for all samples
        3. Running Otter on all samples but one, to infer LIONESS-OTTER networks
    
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

    Notes
    ------

    Toy data:The example gene expression data that we have available here contains gene expression profiles
    for different samples in the columns. Of note, this is just a small subset of a larger gene
    expression dataset. We provided these "toy" data so that the user can test the method.

    References
    ----------
    .. [1]__ 

    Authors: Viola Fanfani
    """

    def __init__(
        self,
        expression_file,
        motif_file,
        ppi_file,
        mode_process="union",
    ):
        """Intialize instance of LionessOtter class and load data."""

        self.expression_file = expression_file
        self.ppi_file = ppi_file
        self.motif_file = motif_file
        self.mode_process = mode_process


        #if not os.path.exists(self.output_folder):
        #    os.makedirs(self.output_folder)

        
        # PPI variables
        self.ppi_data = None
        self.ppi_tfs = None
        
        # gene expression variables
        self.expression_data = None
        self.expression_genes = None
        self.samples = None
        self.n_samples = None

        # motif variables
        self.motif_data = None
        self.motif_tfs = None
        self.motif_genes = None

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
        with Timer("Reading PPI data..."):
            # Read ppi and tfs
            self.ppi_data, self.ppi_tfs = io.read_ppi(self.ppi_file)

        with Timer("Reading expression data..."):
            # Read expression
            self.expression_data, self.expression_genes = io.prepare_expression(
                self.expression_file
            )
            # #TODO: Here fix in case of start/end or list of samples
            self.samples = self.expression_data.columns.tolist()
            self.n_samples = len(self.samples)

        with Timer("Reading Motif data..."):
            # Read the motif
            self.motif_data, self.motif_tfs, self.motif_genes = io.read_motif(
                self.motif_file
            )

        # depending on the strategy for 
        if self.mode_process=='intersection':
            self.universe_genes = sorted(list(set(self.expression_genes).intersection(set(self.motif_genes))))
            self.universe_tfs = sorted(list(set(self.ppi_tfs).intersection(set(self.motif_tfs))))
        else:
            sys.exit('Only intersection is an available modeProcess for the moment')

        # Auxiliary dicts
        self.gene2idx = {x: i for i, x in enumerate(self.universe_genes)}
        self.tf2idx = {x: i for i, x in enumerate(self.universe_tfs)}

        # sort the gene expression, ppi data, motif data
        self.expression_data = self.expression_data.loc[self.universe_genes,self.samples]
        self.ppi_data = self.ppi_data.loc[self.universe_tfs,self.universe_tfs]
        self.motif_data = self.motif_data.loc[self.universe_tfs, self.universe_genes]
        

    def run_otter(self, output_file = 'otter.h5', computing= 'cpu', precision = 'single', lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
        
        """Run LIONESS otter

        Args:
            output_file(str): file to save the otter network. Pass with an extension (txt/h5/csv/npy/mat)
            computing (str, optional): computing for coexpression lioness. Defaults to 'cpu'.
            precision (str, optional): precision. Defaults to single
            lam (float, optional): lambda parameter for otter. Defaults to 0.035
            gamma (float, optional): gamma parameter for otter. Defaults to 0.335
            Iter (int, optional): Iterations for otter. Defaults to 60
            eta (float, optional):  eta parameter for otter. Defaults to 0.00001
            bexp (float, optional): bexp parameter for otter. Defaults to 1
            
        """

        lioness_start = time.time()

        self.save_fmt = output_file.split('.')[-1]
        
        self.output_file = output_file
        
        if precision=='single':
            self.atype = 'float32'
        elif precision=='double':
            self.atype = 'float64'
        else: 
            sys.exit('Precision %s unknonw' %str(precision))
        
        # let's sort the expression and ppi data

        self.expression_data = self.expression_data.astype(self.atype)
        self.motif_data = self.motif_data.astype(self.atype)
        self.ppi_data = self.ppi_data.astype(self.atype)
        
        # First we get the 
        if computing=='cpu':
            correlation_complete = np.corrcoef(self.expression_data.values)
            
            with Timer("Running OTTER for all samples ..."):
                # running otter for all samples
                self.all_otter = otter(self.motif_data.values, self.ppi_data.values, correlation_complete, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
            
        elif computing=='gpu':
            import cupy as cp
            cp._default_memory_pool.free_all_blocks() 
            
            correlation_complete = cp.corrcoef(self.expression_data.values).astype(self.atype)

            with Timer("Running OTTER for all samples ..."):
                # running otter for all samples
                self.all_otter = otter_gpu(self.motif_data.values, self.ppi_data.values, correlation_complete, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
                del correlation_complete
                cp._default_memory_pool.free_all_blocks() 
        else:
            sys.exit('Use one of gpu/cpu for computations')
                
        with Timer("Saving OTTER for all samples ..."):
                   
            print('Saving %s' %self.output_file)
            self.__lioness_to_disk(self.all_otter, path = self.output_file)
                
        print('FINISHED! OTTER took %.2f sec.' % (time.time() - lioness_start))

        
        
    def run_lioness_otter(self, output_folder, save_single = False, save_memory = False, save_fmt = '.h5', online_coexpression = False, computing= 'cpu', cores = 1, precision = 'single', lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
        
        """Run LIONESS otter

        Args:
            output_folder (str): folder to save all the files
            save_memory (bool, optional): whether to save the coexpression with the gene names
            online_coexpression (bool, optional): if true coexpression is computed with a closed form
            computing (str, optional): computing for coexpression lioness. Defaults to 'cpu'.
            cores (int, optional): cores. Defaults to 1.
            precition (str, optional): precision. Defaults to single
            lam (float, optional): lambda parameter for otter. Defaults to 0.035
            gamma (float, optional): gamma parameter for otter. Defaults to 0.335
            Iter (int, optional): Iterations for otter. Defaults to 60
            eta (float, optional):  eta parameter for otter. Defaults to 0.00001
            bexp (float, optional): bexp parameter for otter. Defaults to 1
            
        """

        lioness_start = time.time()

        self.save_fmt = save_fmt
        self.save_single = save_single

        self.output_folder = output_folder
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
        if precision=='single':
            self.atype = 'float32'
        elif precision=='double':
            self.atype = 'float64'
        else: 
            sys.exit('Precision %s unknonw' %str(precision))
        
        # let's sort the expression and ppi data

        self.expression_data = self.expression_data.astype(self.atype)
        self.motif_data = self.motif_data.astype(self.atype)
        self.ppi_data = self.ppi_data.astype(self.atype)
        
        # First we get the 
        if computing=='cpu':
            correlation_complete = np.corrcoef(self.expression_data.values)
            
            with Timer("Running OTTER for all samples ..."):
                # running otter for all samples
                self.all_otter = otter(self.motif_data.values, self.ppi_data.values, correlation_complete, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
            
            with Timer("Saving OTTER for all samples ..."):
                name = "otter.%s" %(str(self.save_fmt))            
                print('Saving %s' %name)
                self.__lioness_to_disk(self.all_otter, path = self.output_folder+'/'+name)
            # Now for each sample we compute the lioness network from correlations and 
            # the panda using the motif and ppi tables
            for s,sample in enumerate(self.samples):
                with Timer("Running OTTER for %d-th sample (%s) ..." %(s,str(sample))):
                    # Run single sample otter with the loop
                    lioness_otter = self._lioness_otter_loop(sample, save_memory=save_memory, computing=computing, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
        
        elif computing=='gpu':
            import cupy as cp
            cp._default_memory_pool.free_all_blocks() 
            
            correlation_complete = cp.corrcoef(self.expression_data.values).astype(self.atype)

            with Timer("Running OTTER for all samples ..."):
                # running otter for all samples
                self.all_otter = otter_gpu(self.motif_data.values, self.ppi_data.values, correlation_complete, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
                del correlation_complete
                cp._default_memory_pool.free_all_blocks() 
                
            with Timer("Saving OTTER for all samples ..."):
                name = "otter.%s" %(str(self.save_fmt))            
                print('Saving %s' %name)
                self.__lioness_to_disk(self.all_otter, path = self.output_folder+'/'+name)
                
            # Now for each sample we compute the lioness network from correlations and 
            # the panda using the motif and ppi tables
            for s,sample in enumerate(self.samples):
                with Timer("Running OTTER for %d-th sample (%s) ..." %(s,str(sample))):
                    # Run single sample otter with the loop
                    lioness_otter = self._lioness_otter_loop(sample, save_memory=save_memory, computing=computing, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)

        print('FINISHED! All LIONESS-OTTER took %.2f sec.' % (time.time() - lioness_start))


    def _lioness_otter_loop(self, sample, save_memory = True, computing = 'cpu', lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
        """Runs LIONESS OTTER for one sample

        Args:
            sample (_type_): _description_
            save_memory (bool, optional): _description_. Defaults to True.
            online_coexpression (bool, optional): _description_. Defaults to False.
            computing (str, optional): _description_. Defaults to 'cpu'.
            computing_panda (str, optional): _description_. Defaults to 'cpu'.
            alpha (float, optional): _description_. Defaults to 0.1.
        """
        
        if not os.path.exists(self.output_folder+'/lioness_otter/'):
            os.makedirs(self.output_folder+'/lioness_otter/')
        
        all_samples_minus_q = list(set(self.samples).difference(set([sample])))
        #correlation_current = self.expression_data.loc[:, all_samples_minus_q ].T.corr()
        
        if computing=='cpu':
            # compute correlation for all samples but q
            correlation_current = np.corrcoef(self.expression_data.loc[:, all_samples_minus_q ].values.astype(self.atype)).astype(self.atype)
            # compute otter for all samples but q
            all_minus_q = otter(self.motif_data.values, self.ppi_data.values, correlation_current, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
        elif computing=='gpu':
            import cupy as cp
            cp._default_memory_pool.free_all_blocks() 
            cp._default_pinned_memory_pool.free_all_blocks() 
            # compute correlation for all samples but q
            correlation_complete = cp.corrcoef(self.expression_data.loc[:, all_samples_minus_q ].values).astype(self.atype)
            # running otter for all samples but q
            all_minus_q = otter_gpu(self.motif_data.values, self.ppi_data.values, correlation_complete, lam=lam, gamma=gamma, Iter=Iter, eta=eta, bexp=bexp)
            del correlation_complete
            
            cp._default_memory_pool.free_all_blocks() 
            cp._default_pinned_memory_pool.free_all_blocks() 
        else:
            sys.exit('Use one of gpu/cpu for computations')
        
        lioness_otter = (self.n_samples*self.all_otter)-(self.n_samples-1)*(all_minus_q)
        
        if self.save_single:
            name = "lioness_otter.%s.%s" %(str(sample), str(self.save_fmt))
            print('Saving %s' %name)
            self.__lioness_to_disk(lioness_otter, path = self.output_folder+'lioness_otter/'+name)
            return(np.array([0]))
        else:
            return(lioness_otter)
        
    
    def _save_single_panda_net(self, net, output_file):

        tab = pd.DataFrame(net, columns = self.universe_genes )
        tab['tf'] = self.universe_tfs

        tab = pd.melt(tab, id_vars='tf', value_vars=tab.columns,var_name='gene', value_name='otter')
        tab['motif'] = self.motif_data.values.flatten(order = 'F')
        tab.to_csv(output_file, sep = '\t', index = False, columns = ['tf', 'gene','motif','otter'])

    def __lioness_to_disk(self, net, path):
        if self.save_fmt == "txt":
            self._save_single_panda_net(net, path)
            np.savetxt(path, net)
        elif self.save_fmt == "h5":
            pd.DataFrame(data = net, columns=self.universe_genes, index = self.universe_tfs).to_hdf(path, key = 'network' ,mode='w')
        elif self.save_fmt == 'csv':
            pd.DataFrame(data = net, columns=self.universe_genes, index = self.universe_tfs).to_csv(path)
        elif self.save_fmt == "npy":
            np.save(path, net)
        elif self.save_fmt == "mat":
            from scipy.io import savemat
            savemat(path, {"PredNet": net})
        else:
            print("Unknown format %s! Use npy format instead." % self.save_fmt)
            np.save(path, net)