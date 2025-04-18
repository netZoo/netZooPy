[![master](https://github.com/netZoo/netZooPy/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/netZoo/netZooPy/actions/workflows/main.yml)
[![devel](https://github.com/netZoo/netZooPy/actions/workflows/main.yml/badge.svg?branch=devel)](https://github.com/netZoo/netZooPy/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/netZoo/netZooPy/branch/devel/graph/badge.svg)](https://codecov.io/gh/netZoo/netZooPy)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/netzoopy/badge/?version=latest)](https://netzoopy.readthedocs.io/en/latest/?badge=latest)
[![tutorials](https://img.shields.io/badge/netZooPy-tutorials-9cf)](https://github.com/netZoo/netZooPy/tree/master/tutorials)
[![alt-version](https://anaconda.org/netzoo/netzoopy/badges/version.svg)](https://anaconda.org/netzoo/netzoopy)
[![Netbooks](https://img.shields.io/badge/netZooPy-netBooks-ff69b4)](http://netbooks.networkmedicine.org/)
[![discussions](https://img.shields.io/badge/netZooPy-discussions-orange)](https://github.com/netZoo/netZooPy/discussions)

netZooPy is tested on: (OS: Ubuntu + Macos) X (Language: Python v3.7 + Python v3.8 + Python v3.9 + Python v3.10)



## Description

netZooPy is a python package to reconstruct, analyse, and plot biological networks.

**WARNING**: for macos arm64 architectures you might have to manually install pytables. We are only testing macos-13
intel architecture for the moment


**WARNING**: the OTTER CLI and class are still relying on a simple approach for reading and merging. Please be careful
if you have NAs and want a non-intersection between W,P,C please rely on PANDA or on your own filtering. 

## Features

netZooPy currently integrates
(gpu)PANDA, (gpu)LIONESS, (gpu)PUMA, SAMBAR, CONDOR, OTTER, DRAGON, COBRA, and BONOBO.

* **PANDA** (Passing Attributes between Networks for Data Assimilation) [[Glass et al. 2013]](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832): PANDA is a method for estimating bipartite gene regulatory networks (GRNs) consisting of two types of nodes: transcription factors (TFs) and genes. An edge between TF $i$ and gene $j$ indicates that gene $j$ is regulated by TF $i$. The edge weight represents the strength of evidence for this regulatory relationship obtained by integrating three types of biological data: gene expression data, protein-protein interaction (PPI) data, and transcription factor binding motif (TFBM) data. PANDA is an iterative approach that begins with a seed GRN estimated from TFBMs and uses message passing between data types to refine the seed network to a final GRN that is consistent with the information contained in gene expression, PPI, and TFBM data. 
  
* **PUMA** (PANDA Using MicroRNA Associations) [[Kuijjer et al.]]("https://www.sciencedirect.com/science/article/pii/S2589004219300872") extends the PANDA framework to model how
microRNAs (miRNAs) participate in gene regulatory networks. PUMA networks are bipartite networks that consist of a
regulatory layer and a layer of genes being regulated, similar to PANDA networks. While the regulatory layer of PANDA
networks consists only of transcription factors (TFs), the regulatory layer of PUMA networks consists of both TFs and
miRNAs. A PUMA network is seeded using a combination of input data sources such as motif scans or ChIP-seq data (for
TF-gene edges) and an miRNA target prediction tool such as TargetScan or miRanda (for miRNA-gene edges). PUMA uses a
message passing framework similar to PANDA to integrate this prior information with gene-gene coexpression and
protein-protein interactions to estimate a final regulatory network incorporating miRNAs. Kuijjer and colleagues [7]
apply PUMA to 38 GTEx tissues and demonstrate that PUMA can identify important patterns in tissue-specific regulation of
genes by miRNA.

* **CONDOR** (COmplex Network Description Of Regulators) [[Platig et al. 2016]](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033): CONDOR is a tool for community detection in bipartite networks. Many community detection methods for unipartite networks are based on the concept of maximizing a modularity metric that compares the weight of edges within communities to the weight of edges between communities, prioritizing community assignments with higher values of the former relative to the latter. CONDOR extends this concept to bipartite networks by optimizing a bipartite version of modularity defined by [[Barber (2007)]](https://pubmed.ncbi.nlm.nih.gov/18233893/). To enable bipartite community detection on large networks such gene regulatory networks, CONDOR uses a fast unipartite modularity maximization method on one of the two unipartite projections of the bipartite network.  In Platig et al. (2016), CONDOR is applied to bipartite networks of single nucleotide polymorphisms (SNPs) and gene expression, where a network edge from a SNP node to a gene node is indicative of an association between the SNP and the gene expression level, commonly known as an expression quantitative trait locus (eQTL). Communities detected with CONDOR contained local hub nodes ("core SNPs") enriched for association with disease, suggesting that functional eQTL relationships are encoded at the community level.

* **LIONESS** (Linear Interpolation to Obtain Network Estimates for Single Samples) [[Kuijjer et al. 2019]](https://doi.org/10.1016/j.isci.2019.03.021): LIONESS is a flexible method for single-sample network integration. The machinery behind LIONESS is a leave-one-out approach. To construct a single-sample network for sample $i$, a first network is estimated on the full dataset and a second network is estimated on the dataset with sample $i$ withheld. The single-sample network is then estimated based on the difference between these two networks. Any method that can be used to estimate a network can be used with LIONESS to estimate single-sample networks. Two common use cases are the use of LIONESS to generate single-sample GRNs based on PANDA and the use of LIONESS to generate single-sample Pearson correlation networks.

* **SAMBAR** (Subtyping Agglomerated Mutations By Annotation Relations) [[Kuijjer et al.]](https://www.nature.com/articles/s41416-018-0109-7): SAMBAR is a tool for studying cancer subtypes based on patterns of somatic mutations in curated biological pathways. Rather than characterize cancer according to mutations at the gene level, SAMBAR agglomerates mutations within pathways to define a pathway mutation score. To avoid bias based on pathway representation, these pathway mutation scores correct for the number of genes in each pathway as well as the number of times each gene is represented in the universe of pathways. By taking a pathway rather than gene-by-gene lens, SAMBAR both de-sparsifies somatic mutation data and incorporates important prior biological knowledge. Kuijjer et al. (2018) demonstrate that SAMBAR is capable of outperforming other methods for cancer subtyping, producing subtypes with greater between-subtype distances; the authors use SAMBAR for a pan-cancer subtyping analysis that identifies four diverse pan-cancer subtypes linked to distinct molecular processes. 

* **OTTER** (Optimization to Estimate Regulation) [[Weighill et al.]](https://www.biorxiv.org/content/10.1101/2020.06.23.167999v2.abstract): OTTER is a GRN inference method based on the idea that observed biological data (PPI data and gene co-expression data) are projections of a bipartite GRN between TFs and genes. Specifically, PPI data represent the projection of the GRN onto the TF-TF space and gene co-expression data represent the projection of the GRN onto the gene-gene space. OTTER reframes the problem of GRN inference as a problem of relaxed graph matching and finds a GRN that has optimal agreement with the observed PPI and coexpression data. The OTTER objective function is tunable in two ways: first, one can prioritize matching the PPI data or the coexpression data more heavily depending on one's confidence in the data source; second, there is a regularization parameter that can be applied to induce sparsity on the estimated GRN. The OTTER objective function can be solved using spectral decomposition techniques and gradient descent; the latter is shown to be closely related to the PANDA message-passing approach (Glass et al. 2013).
  
**WARNING**: the OTTER CLI and class are still relying on a simple approach for reading and merging. Please be careful
if you have NAs and want a non-intersection between W,P,C please rely on PANDA or on your own filtering. 

* **DRAGON** (Determining Regulatory Associations using Graphical models on Omics Networks) [[Shutta et al.]](https://arxiv.org/abs/2104.01690) is a method for estimating multiomic Gaussian graphical models (GGMs, also known as partial correlation networks) that incorporate two different omics data types. DRAGON builds off of the popular covariance shrinkage method of Ledoit and Wolf with an optimization approach that explicitly accounts for the differences in two separate omics "layers" in the shrinkage estimator. The resulting sparse covariance matrix is then inverted to obtain a precision matrix estimate and a corresponding GGM.  Although GGMs assume normally distributed data, DRAGON can be used on any type of continuous data by transforming data to approximate normality prior to network estimation. Currently, DRAGON can be applied to estimate networks with two different types of omics data. Investigators interested in applying DRAGON to more than two types of omics data can consider estimating pairwise networks and "chaining" them together.

* **COBRA** (Co-expression Batch Reduction Adjustment). Batch effects and other covariates are known to induce spurious associations in co-expression networks and confound differential gene expression analyses. These effects are corrected for using various methods prior to downstream analyses such as the inference of co-expression networks and computing differences between them. In differential co-expression analysis, the pairwise joint distribution of genes is considered rather than independently analyzing the distribution of expression levels for each individual gene. Computing co-expression matrices after standard batch correction on gene expression data is not sufficient to account for the possibility of batch-induced changes in the correlation between genes as existing batch correction methods act solely on the marginal distribution of each gene. Consequently, uncorrected, artifactual differential co-expression can skew the correlation structure such that network-based methods that use gene co-expression can produce false, nonbiological associations even using data corrected using standard batch correction. Co-expression Batch Reduction Adjustment (COBRA) addresses this question by computing a batch-corrected gene co-expression matrix based on estimating a conditional covariance matrix. COBRA estimates a reduced set of parameters that express the co-expression matrix as a function of the sample covariates and can be used to control for continuous and categorical covariates. The method is computationally fast and makes use of the inherently modular structure of genomic data to estimate accurate gene regulatory associations and enable functional analysis for high-dimensional genomic data.


* **BONOBO** (Bayesian Optimized Networks Obtained By assimilating Omics data) is a scalable Bayesian model for deriving individual sample-specific co-expression networks by recognizing variations in molecular interactions across individuals. For every sample, BONOBO assumes a Gaussian distribution on the log-transformed centered gene expression and a conjugate prior distribution on the sample-specific co-expression matrix constructed from all other samples in the data. Combining the sample-specific gene expression with the prior distribution, BONOBO yields a closed-form solution for the posterior distribution of the sample-specific co-expression matrices



## Quick guide

Clone the repository into your local disk:

```bash
git clone https://github.com/netZoo/netZooPy.git
```

Then install netZooPy through pip:

```bash
cd netZooPy
pip3 install -e .
```

Upon completion you can load netZooPy in your python code through

```python
import netZooPy
```

## Conda installation

On anaconda.org you will find the conda recipes for all platforms. We recommend using conda environments to keep your analyses self-contained and reproducible.

To install netzoopy through conda:

```bash
conda install -c netzoo -c conda-forge netzoopy
```

## User guide

Please refer to the [documentation](https://netzoopy.readthedocs.io/en/latest/) website for installation instructions and usage.

## License

The software is free and is licensed under the GNU General License v3.0, see the file [LICENSE](LICENSE) for details.

## Feedback/Issues

Please report any issues to the [issues page](https://github.com/netZoo/netZooPy/issues).

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

Latest version: 0.10.9
