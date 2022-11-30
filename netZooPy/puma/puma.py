from __future__ import print_function

import math
import time
import pandas as pd
import numpy as np
from .timer import Timer
from netZooPy.panda.panda import Panda
import netZooPy.panda.calculations as calc
from netZooPy.puma.calculations import compute_puma
import os

class Puma(object):
    """
    Using PUMA to infer gene regulatory network.
        1. Reading in input data (expression data, motif prior, TF PPI data, miR)
        2. Computing coexpression network
        3. Normalizing networks
        4. Running PUMA algorithm
        5. Writing out PUMA network

    Parameters
    ------------

            expression_file : str
                Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.
            motif_file      : str
                Path to file containing the regulation prior as a tab-separated file without a header. This can be a miRNA-Gene predicted network from TargetScan/miRanda.
                However, this can be combined with transcription factor DNA binding motif data in the form of TF-gene-weight(0/1) to estimate gene regulation by TF and miRNA.
                Alternatively, can be a dataframe with the motif network
                If set to none, the gene coexpression matrix is returned as a result network.
            ppi_file        : str
                Path to file containing the TF PPI data, or pandas dataframe.
                This can be provided as 'None' if no TF data is given and PUMA will estimate a miRNA-Gene networks.
            mir_file        : str
                Path to file containing miRNA list or a list. A standard mir_file can be read as:
                >>> with open(mir_file, "r") as f:  
                >>>     miR = f.read().splitlines() 
            computing       : str
                - 'cpu' uses Central Processing Unit (CPU) to run PANDA.
                - 'gpu' use the Graphical Processing Unit (GPU) to run PANDA.
            precision       : str
                'double' computes the regulatory network in double precision (15 decimal digits).
                'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate.
            save_memory     : bool
                True : removes temporary results from memory. The result network is weighted adjacency matrix of size (nTFs, nGenes).
                False: keeps the temporary files in memory. The result network has 4 columns in the form gene - TF - weight in motif prior - PUMA edge.
            save_tmp        : bool
                Save temporary variables.
            remove_missing  : bool
                Removes the gens and TFs that are not present in one of the priors. Works only if modeProcess='legacy'.
            keep_expression_matrix: bool
                Keeps the input expression matrix in the result Puma object.
            modeProcess     : str
                The input data processing mode.
                - 'legacy': refers to the processing mode in netZooPy<=0.5
                - (Default)'union': takes the union of all TFs and genes across priors and fills the missing genes in the priors with zeros.
                - 'intersection': intersects the input genes and TFs across priors and removes the missing TFs/genes.
            alpha           : float
                Learning rate (default: 0.1)
            start           : int
                first sample of expression (default 1)
            end             : int
                last sample of expression (default None)

    Examples
    --------
        Run the PUMA algorithm, leave out motif and PPI data to use Pearson correlation network:
        from netZooPy.puma.puma import Puma
        puma_obj = Puma('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt','../../tests/ToyData/ToyMiRList.txt')

     

    References
    ----------

    ..[1]__ Kuijjer, Marieke L., et al. "PUMA: PANDA Using MicroRNA Associations." OUP Bioinformatics (2019).

        Authors:cychen, davidvi, alessandromarin
    """

    def __init__(
        self,
        expression_file,
        motif_file,
        ppi_file,
        mir_file,
        modeProcess="union",
        computing="cpu",
        precision="double",
        save_memory=False,
        save_tmp=False,
        remove_missing=False,
        keep_expression_matrix=False,
        alpha=0.1,
        start=1,
        end=None,
        df_correlation_matrix = None
    ):
        """
            Intialize instance of Puma class and load data.
        """
        # =====================================================================
        # Data loading
        # =====================================================================
        Panda.processData(
            self,
            modeProcess,
            motif_file,
            expression_file,
            ppi_file,
            remove_missing,
            keep_expression_matrix,
            start=start,
            end=end,
        )

        with Timer("Loading miR data ..."):
            # If the mir_file is a string the mir list is read from file
            # otherwise the input list is used directly
            if type(mir_file) is str:
                with open(mir_file, "r") as f:
                    miR = f.read().splitlines()
            elif isinstance(mir_file,list):
                miR = mir_file
            else:
                raise Exception(
                        "For mir_file please provide either a file name or a list"
                    )
            TFNames = self.unique_tfs
            sort_idx = np.argsort(TFNames)
            self.s1 = sort_idx[np.searchsorted(TFNames, miR, sorter=sort_idx)]

        if remove_missing and motif_file is not None:
            self.__remove_missing()

        # =====================================================================
        # Network construction
        # =====================================================================
        with Timer("Calculating coexpression network ..."):
            if self.expression_data is None:
                if df_correlation_matrix is None:
                    self.correlation_matrix = np.identity(self.num_genes, dtype=int)
                else:
                    #if a gene is missing from the correlation matrix will have NaN in the matrix
                    self.correlation_matrix = df_correlation_matrix.reindex(index=self.gene_names,columns = self.gene_names).values

            else:
                self.correlation_matrix = np.corrcoef(self.expression_data)
            if np.isnan(self.correlation_matrix).any():
                np.fill_diagonal(self.correlation_matrix, 1)
                self.correlation_matrix = np.nan_to_num(self.correlation_matrix)

        if self.motif_data is None:
            print(
                "Returning the correlation matrix of expression data in <Puma_obj>.correlation_matrix"
            )
            # self.puma_network = self.correlation_matrix
            self.__pearson_results_data_frame()
            return
        # Auxiliary dicts
        gene2idx = {x: i for i, x in enumerate(self.gene_names)}
        tf2idx = {x: i for i, x in enumerate(self.unique_tfs)}

        with Timer("Creating motif network ..."):
            self.motif_matrix_unnormalized = np.zeros((self.num_tfs, self.num_genes))
            idx_tfs = [tf2idx[x] for x in self.motif_data[0]]
            idx_genes = [gene2idx[x] for x in self.motif_data[1]]
            idx = np.ravel_multi_index(
                (idx_tfs, idx_genes), self.motif_matrix_unnormalized.shape
            )
            self.motif_matrix_unnormalized.ravel()[idx] = self.motif_data[2]

        if self.ppi_data is None:
            self.ppi_matrix = np.identity(self.num_tfs, dtype=int)
        else:
            with Timer("Creating PPI network ..."):
                self.ppi_matrix = np.identity(self.num_tfs)
                idx_tf1 = [tf2idx[x] for x in self.ppi_data[0]]
                idx_tf2 = [tf2idx[x] for x in self.ppi_data[1]]
                idx = np.ravel_multi_index((idx_tf1, idx_tf2), self.ppi_matrix.shape)
                self.ppi_matrix.ravel()[idx] = self.ppi_data[2]
                idx = np.ravel_multi_index((idx_tf2, idx_tf1), self.ppi_matrix.shape)
                self.ppi_matrix.ravel()[idx] = self.ppi_data[2]

        # =====================================================================
        # Network normalization
        # =====================================================================
        with Timer("Normalizing networks ..."):
            self.correlation_matrix = self._normalize_network(self.correlation_matrix)
            with np.errstate(invalid="ignore"):  # silly warning bothering people
                self.motif_matrix = self._normalize_network(
                    self.motif_matrix_unnormalized
                )
            self.ppi_matrix = self._normalize_network(self.ppi_matrix)
            if precision == "single":
                self.correlation_matrix = np.float32(self.correlation_matrix)
                self.motif_matrix = np.float32(self.motif_matrix)
                self.ppi_matrix = np.float32(self.ppi_matrix)

        # =====================================================================
        # Clean up useless variables to release memory
        # =====================================================================
        if save_memory:
            print("Clearing motif and ppi data, unique tfs, and gene names for speed")
            del (
                self.motif_data,
                self.ppi_data,
                self.unique_tfs,
                self.gene_names,
                self.motif_matrix_unnormalized,
            )

        # =====================================================================
        # Saving middle data to tmp
        # =====================================================================
        if save_tmp:
            with Timer("Saving expression matrix and normalized networks ..."):
                os.makedirs('./tmp',exist_ok=True) 
                if self.expression_data is not None:
                    np.save("./tmp/expression.npy", self.expression_data.values)
                np.save("./tmp/motif.normalized.npy", self.motif_matrix)
                np.save("./tmp/ppi.normalized.npy", self.ppi_matrix)
                np.save("./tmp/correlation_matrix.npy", self.correlation_matrix)

        # Clean up useless variables to release memory
        if keep_expression_matrix:
            self.expression_matrix = self.expression_data.values
        del self.expression_data

        # =====================================================================
        # Running PUMA algorithm
        # =====================================================================
        print("Running PUMA algorithm ...")
        self.puma_network = self.puma_loop(
            self.correlation_matrix,
            self.motif_matrix,
            self.ppi_matrix,
            computing,
            alpha,
        )

    def __remove_missing(self):
        """
            Removes the gens and TFs that are not present in one of the priors. Works only if modeProcess='legacy'.
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
        """
            Standardizes the input data matrices.

        Parameters
        ----------
            x     : array
                Input adjacency matrix.

        Returns
        --------
            normalized_matrix: array
                Standardized adjacency matrix.
        """
        return calc.normalize_network(x)

    def puma_loop(
        self, correlation_matrix, motif_matrix, ppi_matrix, computing="cpu", alpha=0.1
    ):
        """
        The PUMA algorithm.

        Parameters
        ------------
            correlation_matrix: array
                Input coexpression matrix.
            motif_matrix      : array
                Input motif regulation prior network.
            ppi_matrix        : array
                Input PPI matrix.
            computing         : str
                - 'cpu' uses Central Processing Unit (CPU) to run PANDA.
                - 'gpu' use the Graphical Processing Unit (GPU) to run PANDA.

        """

        puma_loop_time = time.time()
        # TODO:This should be using self.correlation. Keeping for retrocompatibility
        motif_matrix = compute_puma(
            correlation_matrix,
            motif_matrix,
            ppi_matrix,
            self.s1,
            computing=computing,
            alpha=alpha,
        )

        print("Running puma took: %.2f seconds!" % (time.time() - puma_loop_time))
        self.puma_network = motif_matrix
        # Ale: reintroducing the export_puma_results array if Puma called with save_memory=False
        if hasattr(self, "unique_tfs"):
            tfs = np.tile(self.unique_tfs, (len(self.gene_names), 1)).flatten()
            genes = np.repeat(self.gene_names, self.num_tfs)
            motif = self.motif_matrix_unnormalized.flatten(order="F")
            force = motif_matrix.flatten(order="F")
            # TODO: should we keep formats consistent between Panda and Puma?
            # self.export_puma_results = pd.DataFrame({'tf':tfs, 'gene': genes,'motif': motif, 'force': force})
            self.export_puma_results = np.column_stack((tfs, genes, motif, force))
        return motif_matrix

    def __pearson_results_data_frame(self):
        """
            Saves PUMA network in edges format.
        """
        genes_1 = np.tile(self.gene_names, (len(self.gene_names), 1)).flatten()
        genes_2 = (
            np.tile(self.gene_names, (len(self.gene_names), 1)).transpose().flatten()
        )
        self.flat_puma_network = self.puma_network.transpose().flatten()
        self.export_puma_results = pd.DataFrame(
            {"tf": genes_1, "gene": genes_2, "force": self.flat_puma_network}
        )
        self.export_puma_results = self.export_puma_results[["tf", "gene", "force"]]
        return None

    def save_puma_results(self, path="puma.npy"):
        """
        Saves PUMA network.

        Parameters
        -------------
            path: str
                Path to save the network.
        """
        with Timer("Saving PUMA network to %s ..." % path):
            # Because there are two modes of operation (save_memory), save to file will be different
            if hasattr(self, "export_puma_results"):
                toexport = self.export_puma_results
            else:
                toexport = self.puma_network
            # Export to file
            if path.endswith(".txt"):
                np.savetxt(path, toexport, fmt="%s", delimiter=" ")
            elif path.endswith(".csv"):
                np.savetxt(path, toexport, fmt="%s", delimiter=",")
            elif path.endswith(".tsv"):
                np.savetxt(path, toexport, fmt="%s", delimiter="\t")
            else:
                np.save(path, toexport)

    def top_network_plot(self, top=100, file="puma_top_100.png"):
        """
        Selects top genes and plot network

        Parameters
        ----------
            top        : int
                Top number of genes to plot.
            file       : str
                File to save the network plot.
        """
        if not hasattr(self, "export_puma_results"):
            raise AttributeError(
                "Puma object does not contain the export_puma_results attribute.\n"
                + "Run Puma with the flag save_memory=False"
            )
        # Ale TODO: work in numpy instead of pandas?
        self.puma_results = pd.DataFrame(
            self.export_puma_results, columns=["tf", "gene", "motif", "force"]
        )
        subset_puma_results = self.puma_results.sort_values(
            by=["force"], ascending=False
        )
        subset_puma_results = subset_puma_results[
            subset_puma_results.tf != subset_puma_results.gene
        ]
        subset_puma_results = subset_puma_results[0:top]
        self.__shape_plot_network(subset_puma_results=subset_puma_results, file=file)
        return None

    def __shape_plot_network(self, subset_puma_results, file="puma.png"):
        """
        Creates plot.

        Paramters:
            subset_puma_results : array
                Reduced PUMA network to the top genes.
            file                : str
                File to save the network plot.
        """
        # reshape data for networkx
        unique_genes = list(
            set(list(subset_puma_results["tf"]) + list(subset_puma_results["gene"]))
        )
        unique_genes = pd.DataFrame(unique_genes)
        unique_genes.columns = ["name"]
        unique_genes["index"] = unique_genes.index
        subset_puma_results = subset_puma_results.merge(
            unique_genes, how="inner", left_on="tf", right_on="name"
        )
        subset_puma_results = subset_puma_results.rename(columns={"index": "tf_index"})
        subset_puma_results = subset_puma_results.drop(["name"], 1)
        subset_puma_results = subset_puma_results.merge(
            unique_genes, how="inner", left_on="gene", right_on="name"
        )
        subset_puma_results = subset_puma_results.rename(
            columns={"index": "gene_index"}
        )
        subset_puma_results = subset_puma_results.drop(["name"], 1)
        links = subset_puma_results[["tf_index", "gene_index", "force"]]
        self.__create_plot(unique_genes=unique_genes, links=links, file=file)
        return None

    def __create_plot(self, unique_genes, links, file="puma.png"):
        """
        Runs the plot.

        Parameters
        -----------
            unique_genes : list
                Unique list of PUMA genes.
            links        : list
                Edges of the subset PUMA network to the top genes.
            file         : str
                File to save the network plot.

        """
        import networkx as nx
        import matplotlib.pyplot as plt

        g = nx.Graph()
        g.clear()
        plt.clf()
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
            """
            Splits the plot label over several lines for plotting purposes.

            Parameters
            ----------
                label: str
                    Input label text.

            Returns
            -------
                label: str
                    Output label text divided over several lines.
            """
            ll = len(label)
            if ll > 6:
                return label[0 : math.ceil(ll / 2)] + "\n" + label[math.ceil(ll / 2) :]
            return label

        for i, l in enumerate(unique_genes.iloc[:, 0]):
            labels[i] = split_label(l)
        pos = nx.spring_layout(g)
        # nx.draw_networkx(g, pos, labels=labels, node_size=40, font_size=3, alpha=0.3, linewidths = 0.5, width =0.5)
        colors = range(len(edges))
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
        nx.draw_spring(g, k=0.25, iterations=50, **options)
        plt.axis("off")
        plt.savefig(file, dpi=300)
        return None

    def return_puma_indegree(self):
        """
        Computes indegree of PUMA network, only if save_memory = False.
        """
        # subset_indegree = self.export_puma_results.loc[:,['gene','force']]
        subset_indegree = self.puma_results.loc[:, ["gene", "force"]]
        self.puma_indegree = subset_indegree.groupby("gene").sum()
        return self.puma_indegree

    def return_puma_outdegree(self):
        """
        Computes outdegree of PUMA network, only if save_memory = False.
        """
        export_puma_results_pd = pd.DataFrame(
            self.export_puma_results, columns=["tf", "gene", "motif", "force"]
        )
        subset_outdegree = export_puma_results_pd.loc[:, ["tf", "force"]]
        self.puma_outdegree = subset_outdegree.groupby("tf").sum()
        return self.puma_outdegree
