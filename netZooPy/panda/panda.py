from __future__ import print_function
import math
import time
import pandas as pd
from .timer import Timer
import numpy as np
from netZooPy.panda import calculations as calc


class Panda(object):
    """ 
    Using PANDA to infer gene regulatory network.
        1. Reading in input data (expression data, motif prior, TF PPI data)
        2. Computing coexpression network
        3. Normalizing networks
        4. Running PANDA algorithm
        5. Writing out PANDA network


    Parameters
    ----------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.
            motif_file : str 
                Path to file containing the transcription factor DNA binding motif data in the form of
                TF-gene-weight(0/1) as a tab-separated file without a header or pandas dataframe.
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
        motif_file,
        ppi_file,
        computing="cpu",
        precision="double",
        save_memory=True,
        save_tmp=True,
        remove_missing=False,
        keep_expression_matrix=False,
        modeProcess="union",
        alpha=0.1,
        start=1,
        end=None,
        with_header=False
    ):
        """ Intialize instance of Panda class and load data.
        """
        # Read data
        
        self.processData(
            modeProcess,
            motif_file,
            expression_file,
            ppi_file,
            remove_missing,
            keep_expression_matrix,
            start=start,
            end=end,
            with_header=with_header
        )
        print(modeProcess,motif_file,expression_file,ppi_file,save_memory,remove_missing,keep_expression_matrix)
        if hasattr(self, "export_panda_results"):
            return

        # =====================================================================
        # Network normalization
        # =====================================================================

        self.precision = precision
        
        with Timer("Normalizing networks ..."):
            self.correlation_matrix = self._normalize_network(self.correlation_matrix)
            with np.errstate(invalid="ignore"):  # silly warning bothering people
                self.motif_matrix = self._normalize_network(
                    self.motif_matrix_unnormalized
                )
            self.ppi_matrix = self._normalize_network(self.ppi_matrix)
            if self.precision == "single":
                self.correlation_matrix = np.float32(self.correlation_matrix)
                self.motif_matrix = np.float32(self.motif_matrix)
                self.ppi_matrix = np.float32(self.ppi_matrix)
        # =====================================================================
        # Clean up useless variables to release memory
        # =====================================================================
        self.tfs, self.genes = self.unique_tfs, self.gene_names

        if (save_memory==True):
            print("Clearing motif and ppi data, unique tfs, and gene names for speed")
            del self.unique_tfs, self.gene_names, self.motif_matrix_unnormalized

        # =====================================================================
        # Saving middle data to tmp
        # =====================================================================
        if save_tmp:
            with Timer("Saving expression matrix and normalized networks ..."):
                if self.expression_data is not None:
                    np.save("/tmp/expression.npy", self.expression_data.values)
                np.save("/tmp/motif.normalized.npy", self.motif_matrix)
                np.save("/tmp/ppi.normalized.npy", self.ppi_matrix)

        # delete expression data
        del self.expression_data

        # =====================================================================
        # Running PANDA algorithm
        # =====================================================================
        if self.motif_data is not None:
            print("Running PANDA algorithm ...")
            self.panda_network = self.panda_loop(
                self.correlation_matrix,
                self.motif_matrix,
                self.ppi_matrix,
                computing,
                alpha,
            )
            # label dataframe
            self.panda_network = pd.DataFrame(
                self.panda_network, index=self.tfs, columns=self.genes
            )
        else:
            self.panda_network = self.correlation_matrix
            self.__pearson_results_data_frame()
            # label dataframe
            self.panda_network = pd.DataFrame(
                self.panda_network, index=self.genes, columns=self.genes
            )

    def __remove_missing(self):
        """ Removes the genes and TFs that are not present in one of the priors. Works only if modeProcess='legacy'.
        """
        if self.expression_data is not None:
            print("Remove expression not in motif:")
            motif_unique_genes = set(self.motif_data[1])
            len_tot = len(self.expression_data)
            self.expression_data = self.expression_data[
                self.expression_data.index.isin(motif_unique_genes)
            ]
            self.gene_names = self.expression_data.index.tolist()
            self.num_genes = len(self.gene_names)
            print(
                "   {} rows removed from the initial {}".format(
                    len_tot - self.num_genes, len_tot
                )
            )
        # if self.motif_data is not None:
        print("Remove motif not in expression data:")
        len_tot = len(self.motif_data)
        self.motif_data = self.motif_data[
            self.motif_data.iloc[:, 1].isin(self.gene_names)
        ]
        self.unique_tfs = sorted(set(self.motif_data[0]))
        self.num_tfs = len(self.unique_tfs)
        print(
            "   {} rows removed from the initial {}".format(
                len_tot - len(self.motif_data), len_tot
            )
        )
        if self.ppi_data is not None:
            print("Remove ppi not in motif:")
            motif_unique_tfs = np.unique(self.motif_data.iloc[:, 0])
            len_tot = len(self.ppi_data)
            self.ppi_data = self.ppi_data[
                self.ppi_data.iloc[:, 0].isin(motif_unique_tfs)
            ]
            self.ppi_data = self.ppi_data[
                self.ppi_data.iloc[:, 1].isin(motif_unique_tfs)
            ]
            print(
                "   {} rows removed from the initial {}".format(
                    len_tot - len(self.ppi_data), len_tot
                )
            )
        return None

    def _normalize_network(self, x):
        """Standardizes the input data matrices.

        Parameters
        ----------
            x     : array
                Input adjacency matrix.

        Returns
        -------
            normalized_matrix: array
                Standardized adjacency matrix.
        """
        normalized_matrix = calc.normalize_network(x)
        return normalized_matrix

    def processData(
        self,
        modeProcess,
        motif_file,
        expression_file,
        ppi_file,
        remove_missing,
        keep_expression_matrix,
        start=1,
        end=None,
        with_header = False,
    ):
        """ Processes data files into data matrices.

        Parameters
        ----------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. 
                By default, the expression file does not have a header, and the cells ares separated by a tab.
                Pass `with_header=True` if the expression data includes the sample names 
            motif_file : str 
                Path to file containing the transcription factor DNA binding motif data in the form of
                TF-gene-weight(0/1) or pandas dataframe.
                If set to none, the gene coexpression matrix is returned as a result network.
            ppi_file : str
                Path to file containing the PPI data. or pandas dataframe. 
                The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix. 
            remove_missing : bool
                Removes the gens and TFs that are not present in one of the priors. Works only if modeProcess='legacy'.
            keep_expression_matrix : bool
                Keeps the input expression matrix in the result Panda object.
            modeProcess : str
                The input data processing mode.
                - 'legacy': refers to the processing mode in netZooPy<=0.5
                - (Default)'union': takes the union of all TFs and genes across priors and fills the missing genes in the priors with zeros.
                - 'intersection': intersects the input genes and TFs across priors and removes the missing TFs/genes.
            with_header: bool
                pass True when the expression file has a header with the sample names
        """

        # if modeProcess=="legacy":
        # =====================================================================
        # Data loading
        # =====================================================================
        if type(motif_file) is str:
            with Timer("Loading motif data ..."):
                self.motif_data = pd.read_csv(motif_file, sep="\t", header=None)
                self.motif_tfs = sorted(set(self.motif_data[0]))
                self.motif_genes = sorted(set(self.motif_data[1]))
                # self.num_tfs = len(self.unique_tfs)
                # print('Unique TFs:', self.num_tfs)
        elif type(motif_file) is not str:
            if motif_file is None:
                self.motif_data = None
                self.motif_genes = []
                self.motif_tfs = []
            else:
                if not isinstance(motif_file, pd.DataFrame):
                    raise Exception(
                        "Please provide a pandas dataframe for motif data with column names as 'source', 'target', and 'weight'."
                    )
                if ("source" not in motif_file.columns) or (
                    "target" not in motif_file.columns
                ):
                    print('renaming motif columns to "source", "target" and "weight" ')
                    motif_file.columns = ["source", "target", "weight"]
                self.motif_data = pd.DataFrame(motif_file.values)
                self.motif_tfs = sorted(set(motif_file["source"]))
                self.motif_genes = sorted(set(motif_file["target"]))
            # self.num_tfs = len(self.unique_tfs)
            # print('Unique TFs:', self.num_tfs)

        if type(expression_file) is str:
            with Timer("Loading expression data ..."):
                if with_header:
                    # Read data with header
                    self.expression_data = pd.read_csv(
                        expression_file, sep="\t", index_col=0
                    )
                else:
                    self.expression_data = pd.read_csv(
                        expression_file, sep="\t", header=None, index_col=0
                    )

                self.expression_data = self.expression_data.iloc[:, (start-1):end]
                self.expression_genes = self.expression_data.index.tolist()
                self.expression_samples = self.expression_data.columns.astype(str)
                # self.num_genes = len(self.gene_names)
                # print('Expression matrix:', self.expression_data.shape)
        elif type(expression_file) is not str:
            if expression_file is not None:
                if not isinstance(expression_file, pd.DataFrame):
                    raise Exception(
                        "Please provide a pandas dataframe for expression data."
                    )
                self.expression_data = expression_file.iloc[:, (start-1):end]  # pd.read_csv(expression_file, sep='\t', header=None, index_col=0)
                self.expression_genes = self.expression_data.index.tolist()
                self.expression_samples = self.expression_data.columns.astype(str)
                # self.num_genes = len(self.gene_names)
                # print('Expression matrix:', self.expression_data.shape)
            else:
                self.gene_names = self.motif_genes
                self.expression_genes = self.motif_genes
                self.num_genes = len(self.gene_names)
                self.expression_data = (
                    None  # pd.DataFrame(np.identity(self.num_genes, dtype=int))
                )
                print(
                    "No Expression data given: correlation matrix will be an identity matrix of size",
                    len(self.motif_genes),
                )

        if len(self.expression_genes) != len(np.unique(self.expression_genes)):
            print(
                "Duplicate gene symbols detected. Consider averaging before running PANDA"
            )

        if type(ppi_file) is str:
            with Timer("Loading PPI data ..."):
                self.ppi_data = pd.read_csv(ppi_file, sep="\t", header=None)
                self.ppi_tfs = sorted(
                    set(pd.concat([self.ppi_data[0], self.ppi_data[1]]))
                )
                print("Number of PPIs:", self.ppi_data.shape[0])
        elif type(ppi_file) is not str:
            if ppi_file is not None:
                if not isinstance(ppi_file, pd.DataFrame):
                    raise Exception("Please provide a pandas dataframe for PPI data.")
                self.ppi_data = ppi_file  # pd.read_csv(ppi_file, sep='\t', header=None)
                self.ppi_tfs = sorted(
                    set(pd.concat([self.ppi_data[0], self.ppi_data[1]]))
                )
                print("Number of PPIs:", self.ppi_data.shape[0])
            else:
                print(
                    "No PPI data given: ppi matrix will be an identity matrix of size",
                    len(self.motif_tfs),
                )
                self.ppi_data = None
                self.ppi_tfs = self.motif_tfs

        if modeProcess == "legacy" and remove_missing and motif_file is not None:
            self.__remove_missing()
            print('new case')
        if modeProcess == "legacy" and remove_missing==False:
            if expression_file is not None:
                self.gene_names = (
                    self.expression_genes
                )  # sorted( np.unique(self.motif_genes +  self.expression_genes ))
            if motif_file is None:
                self.unique_tfs = self.ppi_tfs
            else:
                self.unique_tfs = (
                    self.motif_tfs
                )  # sorted( np.unique(self.ppi_tfs     +  self.motif_tfs ))

        elif modeProcess == "union":
            self.gene_names = sorted(
                np.unique(self.motif_genes + self.expression_genes)
            )
            self.unique_tfs = sorted(np.unique(self.ppi_tfs + self.motif_tfs))

        elif modeProcess == "intersection":
            if motif_file is None:
                self.gene_names = sorted(np.unique(self.expression_genes))
                self.unique_tfs = sorted(np.unique(self.ppi_tfs))
            else:
                self.gene_names = sorted(
                    np.unique(
                        list(
                            set(self.motif_genes).intersection(
                                set(self.expression_genes)
                            )
                        )
                    )
                )
                self.unique_tfs = sorted(
                    np.unique(list(set(self.ppi_tfs).intersection(set(self.motif_tfs))))
                )

        self.num_genes = len(self.gene_names)
        self.num_tfs = len(self.unique_tfs)

        # Auxiliary dicts
        gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}
        if (
            (modeProcess == "union" or modeProcess == "intersection")
            and (self.expression_data is not None)
            and (self.num_genes != 0)
        ):
            # Initialize data & Populate gene expression
            self.expression = np.zeros((self.num_genes, self.expression_data.shape[1]))
            idx_geneEx = [gene2idx.get(x, np.nan) for x in self.expression_genes]
            filtered_genes = [i for (i, v) in zip(self.expression_genes, idx_geneEx) if ~np.isnan(v)]
            idx_geneEx = [x for x in idx_geneEx if str(x) != 'nan']
            self.expression[idx_geneEx, :] = self.expression_data.loc[filtered_genes].values
            self.expression_data = pd.DataFrame(data=self.expression)

        # =====================================================================
        # Network construction
        # =====================================================================
        with Timer("Calculating coexpression network ..."):
            if self.expression_data is None:
                self.correlation_matrix = np.identity(self.num_genes, dtype=int)
            else:
                self.correlation_matrix = np.corrcoef(self.expression_data)
            if np.isnan(self.correlation_matrix).any():
                np.fill_diagonal(self.correlation_matrix, 1)
                self.correlation_matrix = np.nan_to_num(self.correlation_matrix)

        # Clean up useless variables to release memory
        if keep_expression_matrix:
            if self.expression_data is not None:
                self.expression_matrix = self.expression_data.values
            else:
                self.expression_matrix = None

        if self.motif_data is None:
            print(
                "Returning the correlation matrix of expression data in <Panda_obj>.correlation_matrix"
            )
            self.panda_network = self.correlation_matrix
            self.export_panda_results = self.correlation_matrix
            self.motif_matrix = self.motif_data
            self.ppi_matrix = self.ppi_data
            self.__pearson_results_data_frame()
            self.panda_network = pd.DataFrame(
                self.panda_network,
                index=self.expression_genes,
                columns=self.expression_genes,
            )
            return

        with Timer("Creating motif network ..."):
            self.motif_matrix_unnormalized = np.zeros((self.num_tfs, self.num_genes))
            idx_tfs = [tf2idx.get(x, np.nan) for x in self.motif_data[0]]
            idx_genes = [gene2idx.get(x, np.nan) for x in self.motif_data[1]]
            commind1 = ~np.isnan(idx_tfs) & ~np.isnan(idx_genes)
            idx_tfs = [i for (i, v) in zip(idx_tfs, commind1) if v]
            idx_genes = [i for (i, v) in zip(idx_genes, commind1) if v]
            idx = np.ravel_multi_index(
                (idx_tfs, idx_genes), self.motif_matrix_unnormalized.shape
            )
            self.motif_matrix_unnormalized.ravel()[idx] = self.motif_data[2][commind1]

        if self.ppi_data is None:
            self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
        else:
            with Timer("Creating PPI network ..."):
                self.ppi_matrix = np.identity(self.num_tfs)
                idx_tf1 = [tf2idx.get(x, np.nan) for x in self.ppi_data[0]]
                idx_tf2 = [tf2idx.get(x, np.nan) for x in self.ppi_data[1]]
                commind2 = ~np.isnan(idx_tf1) & ~np.isnan(idx_tf2)
                idx_tf1 = [i for (i, v) in zip(idx_tf1, commind2) if v]
                idx_tf2 = [i for (i, v) in zip(idx_tf2, commind2) if v]
                idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
                self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
                self.ppi_matrix.ravel()[idx] = self.ppi_data[2][commind2]

        return

    def panda_loop(
        self, correlation_matrix, motif_matrix, ppi_matrix, computing="cpu", alpha=0.1
    ):
        """ The PANDA algorithm.

        Parameters
        ----------
            correlation_matrix: array
                Input coexpression matrix.
            motif_matrix      : array
                Input motif regulation prior network.
            ppi_matrix        : array
                Input PPI matrix.
            computing         : str
                'cpu' uses Central Processing Unit (CPU) to run PANDA.
                'gpu' use the Graphical Processing Unit (GPU) to run PANDA.
        """


        panda_loop_time = time.time()
        # TODO:This should be using self.correlation. Keeping for retrocompatibility
        motif_matrix = calc.compute_panda(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            computing=computing,
            alpha=alpha,
        )
        print("Running panda took: %.2f seconds!" % (time.time() - panda_loop_time))
        # Ale: reintroducing the export_panda_results array if Panda called with save_memory=False

        if hasattr(self, "unique_tfs"):
            tfs = np.tile(self.unique_tfs, (len(self.gene_names), 1)).flatten()
            genes = np.repeat(self.gene_names, self.num_tfs)
            motif = self.motif_matrix_unnormalized.flatten(order="F")
            force = motif_matrix.flatten(order="F")
            self.export_panda_results = pd.DataFrame(
                {"tf": tfs, "gene": genes, "motif": motif, "force": force}
            )
            # self.export_panda_results = np.column_stack((tfs,genes,motif,force))
        return motif_matrix

    def __pearson_results_data_frame(self):
        """ Saves PANDA network in edges format.
        """
        genes_1 = np.tile(self.gene_names, (len(self.gene_names), 1)).flatten()
        genes_2 = (
            np.tile(self.gene_names, (len(self.gene_names), 1)).transpose().flatten()
        )
        self.flat_panda_network = self.panda_network.transpose().flatten()
        self.export_panda_results = pd.DataFrame(
            {"tf": genes_1, "gene": genes_2, "force": self.flat_panda_network}
        )
        self.export_panda_results = self.export_panda_results[["tf", "gene", "force"]]
        return None

    def save_panda_results(self, path="panda.npy", save_adjacency=False ):
        """ Saves PANDA network.

        Parameters
        ----------
            path: str
                Path to save the network.
            save_adjacency: bool
                if True the output is an adjacency matrix and not the edge list
        """
        with Timer("Saving PANDA network to %s ..." % path):
            # Because there are two modes of operation (save_memory), save to file will be different
            if not hasattr(self, "unique_tfs"):
                toexport = self.panda_network
            else:
                # save the network with names
                toexport = self.export_panda_results
                if save_adjacency:
                    toexport = pd.pivot_table(toexport, values = 'force', index = 'tf', columns='gene', dropna=False)
                    toexport = toexport.reset_index()
            # Export to file
            if path.endswith(".txt"):
                #np.savetxt(path, toexport, fmt="%s", delimiter=" ")
                toexport.to_csv(path, sep=" ", index=False)
            elif path.endswith(".csv"):
                #np.savetxt(path, toexport, fmt="%s", delimiter=",")
                toexport.to_csv(path, sep=",", index=False)
            elif path.endswith(".tsv"):
                #np.savetxt(path, toexport, fmt="%s", delimiter="\t")
                toexport.to_csv(path, sep="\t", index=False)
            else:
                np.save(path, toexport)

    def top_network_plot(self, top=100, file="panda_top_100.png", plot_bipart=False):
        """ Selects top genes.

        Parameters
        ----------
            top        : int
                Top number of genes to plot.
            file       : str
                File to save the network plot.
            plot_bipart: bool
                Plot the network as a bipartite layout.
        """
        if not hasattr(self, "export_panda_results"):
            raise AttributeError(
                "Panda object does not contain the export_panda_results attribute.\n"
                + "Run Panda with the flag save_memory=False"
            )
        # Ale TODO: work in numpy instead of pandas?
        self.panda_results = pd.DataFrame(
            self.export_panda_results, columns=["tf", "gene", "motif", "force"]
        )
        subset_panda_results = self.panda_results.sort_values(
            by=["force"], ascending=False
        )
        subset_panda_results = subset_panda_results[
            subset_panda_results.tf != subset_panda_results.gene
        ]
        subset_panda_results = subset_panda_results[0:top]
        self.__shape_plot_network(
            subset_panda_results=subset_panda_results,
            file=file,
            plot_bipart=plot_bipart,
        )
        return None

    def __shape_plot_network(
        self, subset_panda_results, file="panda.png", plot_bipart=False
    ):
        """ Creates plot.

        Parameters
        -----------
            subset_panda_results : array
                Reduced PANDA network to the top genes.
            file                 : str
                File to save the network plot.
            plot_bipart: bool
                Plot the network as a bipartite layout.
        """
        # reshape data for networkx
        unique_genes = list(
            set(list(subset_panda_results["tf"]) + list(subset_panda_results["gene"]))
        )
        unique_genes = pd.DataFrame(unique_genes)
        unique_genes.columns = ["name"]
        unique_genes["index"] = unique_genes.index
        subset_panda_results = subset_panda_results.merge(
            unique_genes, how="inner", left_on="tf", right_on="name"
        )
        subset_panda_results = subset_panda_results.rename(
            columns={"index": "tf_index"}
        )
        subset_panda_results = subset_panda_results.drop(["name"], 1)
        subset_panda_results = subset_panda_results.merge(
            unique_genes, how="inner", left_on="gene", right_on="name"
        )
        subset_panda_results = subset_panda_results.rename(
            columns={"index": "gene_index"}
        )
        subset_panda_results = subset_panda_results.drop(["name"], 1)
        links = subset_panda_results[["tf_index", "gene_index", "force"]]
        self.__create_plot(
            unique_genes=unique_genes, links=links, file=file, plot_bipart=plot_bipart
        )
        return None

    def __create_plot(self, unique_genes, links, file="panda.png", plot_bipart=False):
        """ Runs the plot.

        Parameters
        ----------
            unique_genes : list
                Unique list of PANDA genes.
            links        : list
                Edges of the subset PANDA network to the top genes.
            file         : str
                File to save the network plot.
            plot_bipart  : bool
                Plot the network as a bipartite layout.

        Notes
        -----
            split_label: Splits the plot label over several lines for plotting purposes.
        """
        import networkx as nx
        import matplotlib.pyplot as plt

        g = nx.Graph()
        g.clear()
        plt.clf()
        # img = plt.imread("../img/panda.jpg")
        # fig, ax = plt.subplots()
        # ax.imshow(img, extent=[0, 400, 0, 300])
        ##ax.plot(x, x, '--', linewidth=5, color='firebrick')
        g.add_nodes_from(unique_genes["index"])
        edges = []
        for i in range(0, len(links)):
            edges = edges + [
                (
                    links.iloc[i]["tf_index"],
                    links.iloc[i]["gene_index"],
                    float(links.iloc[i]["force"]) / 200,
                )
            ]
        g.add_weighted_edges_from(edges)
        labels = {}

        def split_label(label):
            """ Splits the plot label over several lines for plotting purposes.

            Parameters
            ----------
                label: Input label text.

            Returns:
                label: _
                    Output label text divided over several lines.
            """
            ll = len(label)
            if ll > 6:
                return (
                    label[0 : int(np.ceil(ll / 2))]
                    + "\n"
                    + label[int(np.ceil(ll / 2)) :]
                )
            return label

        for i, l in enumerate(unique_genes.iloc[:, 0]):
            labels[i] = split_label(l)
        if not plot_bipart:
            pos = nx.spring_layout(g)
        else:
            pos = nx.drawing.layout.bipartite_layout(g, set(links["tf_index"]))
        # nx.draw_networkx(g, pos, labels=labels, node_size=40, font_size=3, alpha=0.3, linewidth = 0.5, width =0.5)
        print(plot_bipart)
        if not plot_bipart:
            colors = range(len(edges))
        else:
            colors = list(zip(*edges))[-1]

        options = {
            "alpha": 0.7,
            "edge_color": colors,
            "edge_cmap": plt.cm.Blues,
            "node_size": 110,
            "vmin": -100,
            "width": 2,
            "labels": labels,
            "font_weight": "regular",
            "font_size": 3,
            "linewidths": 20,
        }

        nx.draw_networkx(g, pos=pos, **options)
        plt.axis("off")
        plt.savefig(file, dpi=300)
        return None

    def return_panda_indegree(self):
        """ Computes indegree of PANDA network, only if save_memory = False.
        """
        export_panda_results_pd = pd.DataFrame(
            self.export_panda_results, columns=["tf", "gene", "motif", "force"]
        )
        subset_indegree = export_panda_results_pd.loc[:, ["gene", "force"]]
        self.panda_indegree = subset_indegree.groupby("gene").sum()
        return self.panda_indegree

    def return_panda_outdegree(self):
        """ computes outdegree of PANDA network, only if save_memory = False.
        """
        export_panda_results_pd = pd.DataFrame(
            self.export_panda_results, columns=["tf", "gene", "motif", "force"]
        )
        subset_outdegree = export_panda_results_pd.loc[:, ["tf", "force"]]
        self.panda_outdegree = subset_outdegree.groupby("tf").sum()
        return self.panda_outdegree
