.. netZooPy documentation master file, created by
   sphinx-quickstart on Mon Oct  7 21:05:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

Home
===================================

netZooPy is a python package to reconstruct, analyse, and plot biological networks.

Available animals are:

netZooPy currently integrates
-----------------------------------

 (gpu)PANDA, (gpu)LIONESS, (gpu)PUMA, SAMBAR, CONDOR, OTTER, DRAGON, COBRA, and BONOBO.

- **PANDA** (Passing Attributes between Networks for Data Assimilation) `[Glass et al. 2013]
  <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832>`_: PANDA is a method for estimating
  bipartite gene regulatory networks (GRNs) consisting of two types of nodes: transcription factors (TFs) and genes. An
  edge between TF $i$ and gene $j$ indicates that gene $j$ is regulated by TF $i$. The edge weight represents the
  strength of evidence for this regulatory relationship obtained by integrating three types of biological data: gene
  expression data, protein-protein interaction (PPI) data, and transcription factor binding motif (TFBM) data. PANDA is
  an iterative approach that begins with a seed GRN estimated from TFBMs and uses message passing between data types to
  refine the seed network to a final GRN that is consistent with the information contained in gene expression, PPI, and
  TFBM data. 

- **PUMA** (PANDA Using MicroRNA Associations) [[Kuijjer et
  al.]]("https://www.sciencedirect.com/science/article/pii/S2589004219300872") extends the PANDA framework to model how
microRNAs (miRNAs) participate in gene regulatory networks. PUMA networks are bipartite networks that consist of a
regulatory layer and a layer of genes being regulated, similar to PANDA networks. While the regulatory layer of PANDA
networks consists only of transcription factors (TFs), the regulatory layer of PUMA networks consists of both TFs and
miRNAs. A PUMA network is seeded using a combination of input data sources such as motif scans or ChIP-seq data (for
TF-gene edges) and an miRNA target prediction tool such as TargetScan or miRanda (for miRNA-gene edges). PUMA uses a
message passing framework similar to PANDA to integrate this prior information with gene-gene coexpression and
protein-protein interactions to estimate a final regulatory network incorporating miRNAs. Kuijjer and colleagues [7]
apply PUMA to 38 GTEx tissues and demonstrate that PUMA can identify important patterns in tissue-specific regulation of
genes by miRNA.

- **CONDOR** (COmplex Network Description Of Regulators) `[Platig et al. 2016]
  <http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033>`_: CONDOR is a tool for community
  detection in bipartite networks. Many community detection methods for unipartite networks are based on the concept of
  maximizing a modularity metric that compares the weight of edges within communities to the weight of edges between
  communities, prioritizing community assignments with higher values of the former relative to the latter. CONDOR
  extends this concept to bipartite networks by optimizing a bipartite version of modularity defined by [[Barber
  (2007)]](https://pubmed.ncbi.nlm.nih.gov/18233893/). To enable bipartite community detection on large networks such
  gene regulatory networks, CONDOR uses a fast unipartite modularity maximization method on one of the two unipartite
  projections of the bipartite network.  In Platig et al. (2016), CONDOR is applied to bipartite networks of single
  nucleotide polymorphisms (SNPs) and gene expression, where a network edge from a SNP node to a gene node is indicative
  of an association between the SNP and the gene expression level, commonly known as an expression quantitative trait
  locus (eQTL). Communities detected with CONDOR contained local hub nodes ("core SNPs") enriched for association with
  disease, suggesting that functional eQTL relationships are encoded at the community level.

- **LIONESS** (Linear Interpolation to Obtain Network Estimates for Single Samples) `[Kuijjer et al. 2019]
  <https://doi.org/10.1016/j.isci.2019.03.021>`_: LIONESS is a flexible method for single-sample network integration.
  The machinery behind LIONESS is a leave-one-out approach. To construct a single-sample network for sample $i$, a first
  network is estimated on the full dataset and a second network is estimated on the dataset with sample $i$ withheld.
  The single-sample network is then estimated based on the difference between these two networks. Any method that can be
  used to estimate a network can be used with LIONESS to estimate single-sample networks. Two common use cases are the
  use of LIONESS to generate single-sample GRNs based on PANDA and the use of LIONESS to generate single-sample Pearson
  correlation networks.

- **SAMBAR** (Subtyping Agglomerated Mutations By Annotation Relations) `[Kuijjer et al.]
  <https://www.nature.com/articles/s41416-018-0109-7>`_: SAMBAR is a tool for studying cancer subtypes based on patterns
  of somatic mutations in curated biological pathways. Rather than characterize cancer according to mutations at the
  gene level, SAMBAR agglomerates mutations within pathways to define a pathway mutation score. To avoid bias based on
  pathway representation, these pathway mutation scores correct for the number of genes in each pathway as well as the
  number of times each gene is represented in the universe of pathways. By taking a pathway rather than gene-by-gene
  lens, SAMBAR both de-sparsifies somatic mutation data and incorporates important prior biological knowledge. Kuijjer
  et al. (2018) demonstrate that SAMBAR is capable of outperforming other methods for cancer subtyping, producing
  subtypes with greater between-subtype distances; the authors use SAMBAR for a pan-cancer subtyping analysis that
  identifies four diverse pan-cancer subtypes linked to distinct molecular processes. 

- **OTTER** (Optimization to Estimate Regulation) `[Weighill et al.]
  <https://www.biorxiv.org/content/10.1101/2020.06.23.167999v2.abstract>`_: OTTER is a GRN inference method based on the
  idea that observed biological data (PPI data and gene co-expression data) are projections of a bipartite GRN between
  TFs and genes. Specifically, PPI data represent the projection of the GRN onto the TF-TF space and gene co-expression
  data represent the projection of the GRN onto the gene-gene space. OTTER reframes the problem of GRN inference as a
  problem of relaxed graph matching and finds a GRN that has optimal agreement with the observed PPI and coexpression
  data. The OTTER objective function is tunable in two ways: first, one can prioritize matching the PPI data or the
  coexpression data more heavily depending on one's confidence in the data source; second, there is a regularization
  parameter that can be applied to induce sparsity on the estimated GRN. The OTTER objective function can be solved
  using spectral decomposition techniques and gradient descent; the latter is shown to be closely related to the PANDA
  message-passing approach (Glass et al. 2013).

- **BONOBO** (Bayesian Optimized Networks Obtained By assimilating Omics data) is a scalable Bayesian model for deriving
  individual sample-specific co-expression networks by recognizing variations in molecular interactions across
  individuals. For every sample, BONOBO assumes a Gaussian distribution on the log-transformed centered gene expression
  and a conjugate prior distribution on the sample-specific co-expression matrix constructed from all other samples in
  the data. Combining the sample-specific gene expression with the prior distribution, BONOBO yields a closed-form
  solution for the posterior distribution of the sample-specific co-expression matrices.

- **COBRA** (Co-expression Batch Reduction Adjustment). Batch effects and other covariates are known to induce spurious
  associations in co-expression networks and confound differential gene expression analyses. These effects are corrected
  for using various methods prior to downstream analyses such as the inference of co-expression networks and computing
  differences between them. In differential co-expression analysis, the pairwise joint distribution of genes is
  considered rather than independently analyzing the distribution of expression levels for each individual gene.
  Computing co-expression matrices after standard batch correction on gene expression data is not sufficient to account
  for the possibility of batch-induced changes in the correlation between genes as existing batch correction methods act
  solely on the marginal distribution of each gene. Consequently, uncorrected, artifactual differential co-expression
  can skew the correlation structure such that network-based methods that use gene co-expression can produce false,
  nonbiological associations even using data corrected using standard batch correction. Co-expression Batch Reduction
  Adjustment (COBRA) addresses this question by computing a batch-corrected gene co-expression matrix based on
  estimating a conditional covariance matrix. COBRA estimates a reduced set of parameters that express the co-expression
  matrix as a function of the sample covariates and can be used to control for continuous and categorical covariates.
  The method is computationally fast and makes use of the inherently modular structure of genomic data to estimate
  accurate gene regulatory associations and enable functional analysis for high-dimensional genomic data.

Legacy
-------------------------------------

netZooPy has merged of the following repositories and forks:

- `pypanda <https://github.com/davidvi/pypanda>`_ from David Van Ijzendorn, 

- Alessandro Marin's pypanda and Marieke Kuijjer's pypuma `fork <https://github.com/aless80/pypanda>`_, and 

- Cho-Yi Chen's `fork <https://github.com/QuackenbushLab/pypanda>`_.

- Genis Calderer's `pysambar <https://github.com/genis/pysambar>`_ and `pycondor <https://github.com/genis/pycondr>`_


Getting started
=================

.. toctree::

   install/index

.. toctree::

   tutos/index

Documentation
===================

Netzoopy also has a command line interface that allows to run the animal methods from the command line, for example: 

.. code-block::

   netzoopy panda --e expression.txt --m motif.txt --p ppi.txt --o output_panda.txt

Check the documentation below for a full list of options.

.. toctree::
   functions/cli

.. toctree::

   functions/api


.. toctree::

   changelog



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
