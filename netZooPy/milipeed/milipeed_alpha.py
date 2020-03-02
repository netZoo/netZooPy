from __future__ import print_function

import os, os.path,sys
import numpy as np
import pandas as pd
from .timer import Timer
# sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda


class Milipeed(object):
    """
    Description:
       Using MILIPEED to infer single-sample gene regulatory networks.

    Usage:
        1. Reading in input data (expression data, motif prior, methylation data, TF PPI data)
        2. Computing lioness-style coexpression and motif network
        3. Normalizing networks
        4. Running PANDA algorithm
        5. Writing out MILIPEED networks


    Inputs:
       obj: PANDA object, generated with keep_expression_matrix=True.
       obj.motif_matrix: TF DNA motif binding data in tf-by-gene format.

    Authors: 
       dcolinmorgan
    """

    def __init__(self, expression_file, motif_file,methylation_file, ppi_file, start=1, end=None,save_dir='milipeed_output', save_fmt='txt'):
        # =====================================================================
        # Data loading
        # =====================================================================
        if methylation_file is not None and motif_file is not None:
            with Timer('Loading methylation data ...'):
                tmp = pd.read_csv(methylation_file, sep='\t', header=0,index_col=0)
                self.methylation_map = pd.read_csv('~/netZooPy/tests/ToyData/MotifPrior_CGmap.txt', sep='\t', header=None,index_col=2)
                self.mdata=self.methylation_map.merge(tmp,left_index=True, right_index=True)
                self.methylation_subjects = sorted(set(tmp.columns))
                self.methylation_genes = self.mdata[1].tolist()
                self.methylation_tfs = self.mdata[0].tolist()
                self.methylation_data=self.mdata#.drop([0,1],1)
                # print('Methylation matrix:', self.methylation_data.shape)

            with Timer('Loading motif data ...'):
                self.motif_data = pd.read_csv(motif_file, sep='\t', header=None)
                self.motif_tfs = sorted(set(self.motif_data[0]))
                self.motif_genes = sorted(set(self.motif_data[1]))
        else:
            self.methyl_data = None
            print('No Methylation informed motif information given: motif will be uniform and output lioness networks')

        if expression_file:
            with Timer('Loading expression data ...'):
                self.expression_data = pd.read_csv(expression_file, sep='\t', header=0, index_col=0)
                self.expression_genes = self.expression_data.index.tolist()
                self.expression_subjects = self.expression_data.columns
                print('Expression matrix:', self.expression_data.shape)
        else:
            self.gene_names = list(set(self.motif_data[1]))
            self.num_genes = len(self.gene_names)
            self.expression_data = None  # pd.DataFrame(np.identity(self.num_genes, dtype=int))
            print('No Expression data given: correlation matrix will be an identity matrix of size', self.num_genes)

        if ppi_file:
            with Timer('Loading PPI data ...'):
                self.ppi_data = pd.read_csv(ppi_file, sep='\t', header=None)
                self.ppi_tfs = sorted(set(self.ppi_data[0]))
                print('Number of PPIs:', self.ppi_data.shape[0])
        else:
            print('No PPI data given: ppi matrix will be an identity matrix of size', self.num_tfs)
            self.ppi_data = None

        self.subjects   = sorted(np.unique( list(set(self.expression_subjects).intersection(set(self.methylation_subjects)) )))
        self.me_genes   = sorted(np.unique( list(set(self.expression_genes).intersection(set(self.methylation_genes))) ))
        self.gene_names = sorted(np.unique( list(set(self.me_genes).intersection(set(self.motif_genes))) ))
        self.me_tfs     = sorted(np.unique( list(set(self.ppi_tfs).intersection(set(self.methylation_tfs)) )))
        self.unique_tfs = sorted(np.unique( list(set(self.me_tfs).intersection(set(self.motif_tfs)) )))
        self.num_genes  = len(self.gene_names)
        self.num_tfs    = len(self.unique_tfs)
        self.num_subj   = len(self.subjects)
        self.indexes    = range(self.num_subj)[start-1:end] 

        # Initialize expression data
        self.expression = np.zeros((self.num_genes, self.expression_data.shape[1]))
        # Auxiliary dicts
        self.gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        self.subj2idx = {x: i for i, x in enumerate(self.subjects)}
        self.tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}
        # Populate gene expression
        idx_geneEx = [self.gene2idx.get(x, 0) for x in self.expression_genes]
        # self.EEG=pd.DataFrame(dict.items(gene2idx))
        idx_subjEx = [self.subj2idx.get(x, 0) for x in self.expression_subjects]
        # self.EES=pd.DataFrame(dict.items(subj2idx))

        print(self.expression.shape)
        print(self.expression_data.shape)
        self.expression[idx_geneEx, :] = self.expression_data.values
        self.expression_data = pd.DataFrame(data=self.expression)

        # Initialize methylation data
        self.methylation = np.zeros((self.num_genes, self.methylation_data.shape[1]))
        # Auxiliary dicts
        # gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        # tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}
        # Populate methylation
        idx_geneMe = [self.gene2idx.get(x, 0) for x in self.methylation_genes]
        self.MEG=pd.DataFrame(dict.items(self.gene2idx))[0]
        idx_subjMe = [self.subj2idx.get(x, 0) for x in self.methylation_subjects]
        self.MES=pd.DataFrame(dict.items(self.subj2idx))[0]
        idx_tfMe = [self.tf2idx.get(x, 0) for x in self.methylation_tfs]
        self.MET=pd.DataFrame(dict.items(self.tf2idx))[0]
        print(self.methylation.shape)
        print(self.methylation_data.shape)
        self.methylation[idx_geneMe, :] = self.methylation_data.values
        self.methylation_data = pd.DataFrame(data=self.methylation)
    
        # Create the output folder if not exists
        self.save_dir = save_dir
        self.save_fmt = save_fmt
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Run MILIPEED
        self.total_lioness_network = self.__milipeed_loop()

        # # create result data frame
        self.export_milipeed_results = pd.DataFrame(total_lioness_network)

    def __milipeed_loop(self):
        for iii in self.indexes:
            print("Running MILIPEED for subject %d:" % (iii+1))
            idx = [x for x in range(self.num_subj) if x != iii]  # all samples except iii
            # =====================================================================
            # Network construction
            # =====================================================================
            with Timer('Calculating coexpression network ...'):
                if self.expression_data is None:
                    self.correlation_matrix = np.identity(self.num_genes, dtype=int)
                else:
                    self.correlation_matrix = np.corrcoef(self.expression_data)
                if np.isnan(self.correlation_matrix).any():
                    np.fill_diagonal(self.correlation_matrix, 1)
                    self.correlation_matrix = np.nan_to_num(self.correlation_matrix)
                # correlation_matrix = self._normalize_network(correlation_matrix)
                # self.subset_correlation_network = self._normalize_network(np.corrcoef(self.expression_data.values[:, idx]))
                self.subset_correlation_network = np.corrcoef(self.expression_data.values[:, idx])
                self.lionish_network = self.num_subj * (self.correlation_matrix - self.subset_correlation_network) + self.subset_correlation_network
            # return self
            with Timer('Creating motif network ...'):
                self.motif_matrix_unnormalized = np.zeros((self.num_tfs, self.num_genes))
                idx_tfs = [self.tf2idx.get(x, 0) for x in self.MET]
                idx_genes = [self.gene2idx.get(x, 0) for x in self.MEG]
                idx = np.ravel_multi_index((idx_tfs, idx_genes), self.motif_matrix_unnormalized.shape)
                if self.methylation_data is not None:
                    self.motif_matrix_unnormalized.ravel()[idx] = 1-(self.methylation_data[iii])/100
                else:
                    self.motif_matrix_unnormalized.ravel()[idx] = self.motif_data[2]
            
            if self.ppi_data is None:
                self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
            else:
                with Timer('Creating PPI network ...'):
                    self.ppi_matrix = np.identity(self.num_tfs)
                    idx_tf1 = [self.tf2idx.get(x, 0) for x in self.ppi_data[0]]
                    idx_tf2 = [self.tf2idx.get(x, 0) for x in self.ppi_data[1]]
                    idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                    idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
                    self.ppi_matrix.ravel()[idx] = self.ppi_data[2]

            milipeed_network = self.panda_loop(self.lionish_network, self.motif_matrix_unnormalized, self.ppi_matrix)

            with Timer("Saving MILIPEED network %d to %s using %s format:" % (iii+1, self.save_dir, self.save_fmt)):
                path = os.path.join(self.save_dir, "milipeed.%d.%s" % (iii+1, self.save_fmt))
                if self.save_fmt == 'txt':
                    np.savetxt(path, milipeed_network)
                elif self.save_fmt == 'npy':
                    np.save(path, milipeed_network)
                elif self.save_fmt == 'mat':
                    from scipy.io import savemat
                    savemat(path, {'PredNet': milipeed_network})
                else:
                    print("Unknown format %s! Use npy format instead." % self.save_fmt)
                    np.save(path, milipeed_network)
            if iii == 0:
                self.total_milipeed_network = np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)
            else:
                self.total_milipeed_network=np.column_stack((self.total_milipeed_network ,np.fromstring(np.transpose(milipeed_network).tostring(),dtype=milipeed_network.dtype)))

        return self.total_milipeed_network

    def save_milipeed_results(self, file='milipeed.txt'):
        '''Write milipeed results to file.'''
        np.savetxt(file, self.total_milipeed_network, delimiter="\t",header="")
        return None
