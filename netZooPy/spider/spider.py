from __future__ import print_function

import time

import numpy as np
import pandas as pd

from netZooPy.panda.panda import Panda
import netZooPy.panda.calculations as calc
from netZooPy.panda.timer import Timer


class Spider(Panda):
    """
    Using SPIDER to infer epigenetically-informed gene regulatory networks.

    SPIDER (Seeding PANDA Interactions to Derive Epigenetic Regulation) extends
    PANDA by seeding the message-passing algorithm with a motif prior that has
    been filtered by an epigenetic mask (e.g. open-chromatin / ATAC-seq or
    DNase-seq). Only motif interactions that fall in accessible chromatin are
    retained; the rest are removed before running PANDA's message passing. A
    degree-adjustment step reduces the penalty on hub nodes during the z-score
    normalization.

    Steps:
        1. Read input data (expression, motif prior, epigenetic filter, PPI).
        2. Apply the epigenetic filter as a mask on the motif prior.
        3. Compute the coexpression network.
        4. Degree-adjust and normalize the networks.
        5. Run the PANDA message-passing algorithm.

    Parameters
    ------------
        expression_file : str
            Path to file containing the gene expression data or a pandas
            DataFrame (genes x samples). If None, an identity coexpression
            matrix is used.
        motif_file : str
            Path to the motif prior (tab-separated, no header): TF, gene, weight.
            May also be a pandas DataFrame.
        epifilter_file : str
            Path to a binary epigenetic filter of the SAME rows/order as the
            motif prior (tab-separated, no header): TF, gene, {0,1}. May also be
            a pandas DataFrame. Motif edges whose filter value is 0 are removed.
            If None, SPIDER reduces exactly to PANDA.
        ppi_file : str
            Path to the TF-TF PPI prior (tab-separated, no header): TF1, TF2,
            weight. May also be a pandas DataFrame. If None, the identity is used.
        computing : str
            'cpu' (default) or 'gpu'.
        precision : str
            'double' (default) or 'single'.
        save_memory : bool
            If True, removes intermediate matrices; the result is a weighted
            (nTFs x nGenes) adjacency. If False (default), keeps the edge-list
            export in `export_spider_results`.
        save_tmp : bool
            Save intermediate matrices to a tmp folder.
        remove_missing : bool
            Remove genes/TFs missing from a prior (only for modeProcess='legacy').
        keep_expression_matrix : bool
            Keep the input expression matrix on the resulting object (needed for
            downstream LIONESS).
        modeProcess : str
            'union' (default), 'intersection' or 'legacy'.
        alpha : float
            Learning rate (default 0.1).

    Attributes
    -----------
        spider_network : np.ndarray
            The inferred regulatory network (nTFs x nGenes).
        export_spider_results : np.ndarray or pd.DataFrame
            Long-format edges (only when save_memory=False).

    Examples
    --------
        >>> from netZooPy.spider.spider import Spider
        >>> spider_obj = Spider('expression.txt', 'motif.txt',
        ...                     'epifilter.txt', 'ppi.txt', save_memory=False)
        >>> spider_obj.save_spider_results('spider.txt')

    References
    ----------
    .. [1] Sonawane, Abhijeet Rajendra, et al. "Constructing gene regulatory
       networks using epigenetic data." npj Systems Biology and Applications
       7.1 (2021): 1-13.
    """

    def __init__(
        self,
        expression_file,
        motif_file,
        epifilter_file,
        ppi_file,
        computing="cpu",
        precision="double",
        save_memory=False,
        save_tmp=False,
        remove_missing=False,
        keep_expression_matrix=False,
        modeProcess="union",
        alpha=0.1,
    ):
        # ------------------------------------------------------------------
        # 1. Apply the epigenetic filter to the motif prior BEFORE PANDA reads
        #    it. The filter must match the motif rows one-to-one.
        # ------------------------------------------------------------------
        motif_file = self._apply_epifilter(motif_file, epifilter_file)

        # ------------------------------------------------------------------
        # 2. Load and align data via PANDA's shared processing routine.
        # ------------------------------------------------------------------
        Panda.processData(
            self,
            modeProcess,
            motif_file,
            expression_file,
            ppi_file,
            remove_missing,
            keep_expression_matrix,
        )

        if self.motif_data is None:
            raise ValueError(
                "SPIDER requires a motif prior (with an epigenetic filter)."
            )

        # ------------------------------------------------------------------
        # 3. Coexpression network.
        # ------------------------------------------------------------------
        with Timer("Calculating coexpression network ..."):
            if self.expression_data is None:
                self.correlation_matrix = np.identity(self.num_genes, dtype=int)
            else:
                self.correlation_matrix = np.corrcoef(self.expression_data)
            if np.isnan(self.correlation_matrix).any():
                np.fill_diagonal(self.correlation_matrix, 1)
                self.correlation_matrix = np.nan_to_num(self.correlation_matrix)

        # ------------------------------------------------------------------
        # 4. Build the motif and PPI matrices.
        # ------------------------------------------------------------------
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

        # ------------------------------------------------------------------
        # 5. SPIDER-specific degree adjustment, then normalize.
        # ------------------------------------------------------------------
        with Timer("Degree-adjusting and normalizing networks ..."):
            self.motif_matrix_unnormalized = self._degree_adjust(
                self.motif_matrix_unnormalized
            )
            self.correlation_matrix = calc.normalize_network(self.correlation_matrix)
            with np.errstate(invalid="ignore"):
                self.motif_matrix = calc.normalize_network(
                    self.motif_matrix_unnormalized
                )
            self.ppi_matrix = calc.normalize_network(self.ppi_matrix)
            if precision == "single":
                self.correlation_matrix = np.float32(self.correlation_matrix)
                self.motif_matrix = np.float32(self.motif_matrix)
                self.ppi_matrix = np.float32(self.ppi_matrix)

        if save_memory:
            print("Clearing motif and ppi data, unique tfs, and gene names for speed")
            del (
                self.motif_data,
                self.ppi_data,
                self.unique_tfs,
                self.gene_names,
                self.motif_matrix_unnormalized,
            )

        if save_tmp:
            import os

            with Timer("Saving expression matrix and normalized networks ..."):
                os.makedirs("./tmp", exist_ok=True)
                if self.expression_data is not None:
                    np.save("./tmp/expression.npy", self.expression_data.values)
                np.save("./tmp/motif.normalized.npy", self.motif_matrix)
                np.save("./tmp/ppi.normalized.npy", self.ppi_matrix)

        if keep_expression_matrix:
            self.expression_matrix = self.expression_data.values
        del self.expression_data

        # ------------------------------------------------------------------
        # 6. Run the PANDA message-passing algorithm.
        # ------------------------------------------------------------------
        print("Running SPIDER algorithm ...")
        self.spider_network = self.spider_loop(
            self.correlation_matrix,
            self.motif_matrix,
            self.ppi_matrix,
            computing=computing,
            alpha=alpha,
        )
        # keep PANDA-compatible attribute name too
        self.panda_network = self.spider_network

    # ---------------------------------------------------------------------
    # SPIDER-specific helpers
    # ---------------------------------------------------------------------
    @staticmethod
    def _read_prior(x):
        """Read a prior given as a file path or return a copy of a DataFrame."""
        if isinstance(x, str):
            return pd.read_csv(x, sep="\t", header=None)
        return x.copy()

    @classmethod
    def _apply_epifilter(cls, motif_file, epifilter_file):
        """
        Multiply the motif weight (column 3) by the epigenetic filter (column 3).

        The filter must have the same rows, in the same order, as the motif
        (columns 1 and 2 must match). If epifilter_file is None, the motif is
        returned unchanged (SPIDER == PANDA).
        """
        if epifilter_file is None:
            return motif_file
        motif = cls._read_prior(motif_file)
        epi = cls._read_prior(epifilter_file)
        motif = motif.reset_index(drop=True)
        epi = epi.reset_index(drop=True)
        if motif.shape[0] != epi.shape[0]:
            raise ValueError(
                "Chromatin accessibility data does not match motif data size."
            )
        # columns 0 and 1 (TF, gene) must line up
        if not (
            motif.iloc[:, 0].reset_index(drop=True).equals(
                epi.iloc[:, 0].reset_index(drop=True)
            )
            and motif.iloc[:, 1].reset_index(drop=True).equals(
                epi.iloc[:, 1].reset_index(drop=True)
            )
        ):
            raise ValueError(
                "Chromatin accessibility data does not match motif data order."
            )
        filtered = motif.copy()
        filtered.iloc[:, 2] = motif.iloc[:, 2].to_numpy() * epi.iloc[:, 2].to_numpy()
        return filtered

    @staticmethod
    def _degree_adjust(A):
        """
        Degree adjustment so hub nodes are not penalized by z-score scaling.

        Mirrors netZooR::degreeAdjust:
            k1 = colSums(A)/nrow(A);  k2 = rowSums(A)/ncol(A)
            B  = (k1 broadcast over rows)^2 + (k2 broadcast over cols)^2
            A  = A * sqrt(B)
        """
        A = np.asarray(A, dtype=float)
        n_rows, n_cols = A.shape
        k1 = A.sum(axis=0) / n_rows          # length n_cols
        k2 = A.sum(axis=1) / n_cols          # length n_rows
        B = np.tile(k1, (n_rows, 1)) ** 2
        B = B + np.tile(k2.reshape(-1, 1), (1, n_cols)) ** 2
        return A * np.sqrt(B)

    def spider_loop(
        self, correlation_matrix, motif_matrix, ppi_matrix, computing="cpu", alpha=0.1
    ):
        """Run the shared PANDA message-passing loop on the seeded networks."""
        t0 = time.time()
        network = calc.compute_panda(
            correlation_matrix,
            ppi_matrix,
            motif_matrix,
            computing=computing,
            alpha=alpha,
        )
        print("Running SPIDER took: %.2f seconds!" % (time.time() - t0))
        # Build the long-format edge export when data was retained.
        if hasattr(self, "unique_tfs"):
            tfs = np.tile(self.unique_tfs, (len(self.gene_names), 1)).flatten()
            genes = np.repeat(self.gene_names, self.num_tfs)
            motif = self.motif_matrix_unnormalized.flatten(order="F")
            force = network.flatten(order="F")
            self.export_spider_results = np.column_stack((tfs, genes, motif, force))
        return network

    def save_spider_results(self, path="spider.npy"):
        """Save the SPIDER network. Format is inferred from the extension."""
        with Timer("Saving SPIDER network to %s ..." % path):
            if hasattr(self, "export_spider_results"):
                to_export = self.export_spider_results
            else:
                to_export = self.spider_network
            if path.endswith(".txt"):
                np.savetxt(path, to_export, fmt="%s", delimiter=" ")
            elif path.endswith(".csv"):
                np.savetxt(path, to_export, fmt="%s", delimiter=",")
            elif path.endswith(".tsv"):
                np.savetxt(path, to_export, fmt="%s", delimiter="\t")
            else:
                np.save(path, to_export)
