## Description
This repo is based on the following repos:
- [https://github.com/aless80/pypanda](https://github.com/aless80/pypanda),
- [https://github.com/QuackenbushLab/pypanda](https://github.com/QuackenbushLab/pypanda), 
- which was based on [https://github.com/davidvi/pypanda](https://github.com/davidvi/pypanda).  
- Compared to QuackenbushLab/pypanda this repository adds the Python implementation of PUMA ([run_puma.py](netZooPy/netZooPy/pypuma/run_puma.py) and [netZooPy/netZooPy/pypuma/puma.py](pypanda/puma.py)). 
NaN values in normalized matrices are replaced with values normalized by the overall z-score. This allows running the Toy Data provided in this repository.   
  
## Table of Contents
* [Links to literature](#links-to-literature)
* [Panda algorithm](#panda-algorithm)  
* [Installation](#installation)  
* [Usage](#usage)  
  * [Run from terminal](#run-from-terminal)
  * [Run from python](#run-from-python)
* [Toy data](#toy-data)
* [Results](#results)


## Links to literature 

* **[PUMA](https://static-content.springer.com/esm/art%3A10.1186%2Fs13045-017-0465-4/MediaObjects/13045_2017_465_MOESM3_ESM.pdf)** (PANDA Using MicroRNA Associations)  
_Manuscript in preparation, used in [PUMA](https://static-content.springer.com/esm/art%3A10.1186%2Fs13045-017-0465-4/MediaObjects/13045_2017_465_MOESM3_ESM.pdf)._  
C and MATLAB code: [https://github.com/mararie/PUMA](https://github.com/mararie/PUMA)

* **[PANDA](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832)** Passing Attributes between Networks for Data Assimilation  
_Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions, PLoS One, 2013 May 31;8(5):e64832_  
Original PANDA C++ code: [http://sourceforge.net/projects/panda-net/](http://sourceforge.net/projects/panda-net/).  

* **[LIONESS](https://arxiv.org/abs/1505.06440)** (Linear Interpolation to Obtain Network Estimates for Single Samples)
* **[LIONESSR](https://doi.org/10.1016/j.isci.2019.03.02)** The R version   
_Marieke Lydia Kuijjer, Matthew Tung,GuoCheng Yuan,John Quackenbush, Kimberly Glass. Estimating sample-specific regulatory networks_  

LIONESS can be used to estimate single-sample networks using aggregate networks made with any network reconstruction algorithm (http://arxiv.org/pdf/1505.06440.pdf).

## Panda algorithm
<img src="img/panda.png" height="300">  

To find agreement between the three input networks first the responsibility (R) is calculated.  

<img src="img/responsibility.png" height="30">  

Thereafter availability (A) is calculated.  

<img src="img/availability.png" height="30">  

Availability and responsibility are combined with the following formula.  

<img src="img/combine.png" height="30">  

Protein cooperativity and gene co-regulatory networks are updated.  

<img src="img/cooperativity.png" height="30">  
<img src="/img/co-regulatory.png" height="30">  

P and C are updated to satisfy convergence.  

<img src="img/p.png" height="30">  
<img src="/img/c.png" height="30">  

Hamming distance is calculated every iteration.  

<img src="img/hamming.png" height="40">  


## Installation
PyPanda runs on Python 3. You can either run the pypanda script directly (see [Usage](#usage)) or install it on your system. We recommend the following commands to install pypandas on UNIX systems:
#### Using  a virtual environment
Using [python virtual environments](http://docs.python-guide.org/en/latest/dev/virtualenvs/) is the cleanest installation method. 

Cloning git and setting up a [python virtual environment](http://docs.python-guide.org/en/latest/dev/virtualenvs/):
```no-highlight
pip install --user pipenv   #Make sure you have pipenv
git clone https://github.com/netZoo/netZooPy.git
cd netZooPy
```
Creating a virtual environment and installing pypanda:
```no-highlight
virtualenv pypandaenv #virtual environment created in a folder inside the git folder 
source pypandaenv/bin/activate
(pypandaenv)$ pip3 install -r requirements.txt
(pypandaenv)$ python3 setup.py install --record files.txt
```
Uninstall pypanda from virtual environment:
```no-highlight
cat files.txt | xargs rm -rf
```
Complete removal of virtual environment and pypanda:
```no-highlight
(pypanda)$ deactivate	#Quit virtual environment
rm -rf pypandaenv
```

#### Using pip 
Never use ~~sudo pip~~. Instead you can use pip on the user's install directory:
```no-highlight
git clone https://github.com/netZooPy/netZooPy.git
cd pypanda
python3 setup.py install --user
#to run from the command line you will need to make pypanda executable and add the bin directory to your PATH:
cd bin
chmod +x pypanda
echo "$(pwd):PATH" >> ~/.bashrc
source ~/.bashrc
```

```no-highlight
pip uninstall pypanda
```
To run pypanda from Windows (tested on Windows 10) install Git (https://git-scm.com/downloads) and Anaconda Python3 (https://www.continuum.io/downloads) and from the Anaconda prompt run:
```no-highlight
git clone https://github.com/netZooPy/netZooPy.git
cd netZooPy
python3 setup.py install
```

## Usage
#### Run from terminal
pypanda (or pypuma) can be run directly from the terminal with the following options:
```
-h help
-e, --expression: expression values
-m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided 
-p, --ppi: pair file of PPI edges
-o, --output: output file
-i, --mir: mir data miR file (only for pypuma)
-r, --rm_missing
-q, --lioness: output for Lioness single sample networks 
```
To run pypanda on toy data:
```
python run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o output_panda.txt
```
To reconstruct a single sample Lioness Pearson correlation network (this can take some time):
```python
python3 run_panda.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -o output_panda.txt -q output_lioness.txt
```

To run pypuma on toy data:
```python
python3 run_puma.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -o output_puma.txt -i ../../tests//ToyData/ToyMiRList.txt
```
To reconstruct a single sample Lioness Pearson correlation network using pypuma (this can take some time):
```python
python3 run_puma.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -i ../../tests/ToyData/ToyMiRList.txt -o output_puma.txt -q output_lioness.txt
```
For pypuma see also [PyPuma](https://github.com/netZooPy/netZooPy/pypuma#installation). 

#### Run from python
Fire up your python shell or ipython notebook. Use the python installation in the virtual environment if you installed pypanda there. 

Import the classes in the pypanda library:
```python
from netZooPy.panda.panda import Panda
from netZooPy.puma.puma import Puma
from netZooPy.lioness.lioness import Lioness
```
Run the Panda algorithm, leave out motif and PPI data to use Pearson correlation network:
```python
panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False)
```
Save the results:
```python
panda_obj.save_panda_results('Toy_Panda.pairs.txt')
```
Return a network plot:

```python
panda_obj.top_network_plot(top=70, file='top_genes.png')
```
Calculate in- and outdegrees for further analysis:
```python
indegree = panda_obj.return_panda_indegree()
outdegree = panda_obj.return_panda_outdegree()
```
To run the Lioness algorithm for single sample networks, first run panda (or puma) using the keep_expression_matrix flag, then use Lioness as follows:
```python
panda_obj = Panda('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt', remove_missing=False, keep_expression_matrix=True)
lioness_obj = Lioness(panda_obj)
```
Save Lioness results:
```python
lioness_obj.save_lioness_results('Toy_Lioness.txt')
```
Return a network plot for one of the Lioness single sample networks:
```python
plot = AnalyzeLioness(lioness_obj)
plot.top_network_plot(column= 0, top=100, file='top_100_genes.png')
```

Run the Puma algorithm, leave out motif and PPI data to use Pearson correlation network:
```python
puma_obj = Puma('../../tests/ToyData/ToyExpressionData.txt', '../../tests/ToyData/ToyMotifData.txt', '../../tests/ToyData/ToyPPIData.txt','../../tests/ToyData/ToyMiRList.txt')
```

## Toy data
The example gene expression data that we have available here contains gene expression profiles for different samples in the columns. Of note, this is just a small subset of a larger gene expression dataset. We provided these "toy" data so that the user can test the method. 

However, if you plan to model gene regulatory networks on your own dataset, you should use your own expression data as input.

## Results
```
Example Panda output:
TF  Gene  Motif Force
---------------------
CEBPA	AACSL	0.0	-0.951416589143
CREB1	AACSL	0.0	-0.904241609324
DDIT3	AACSL	0.0	-0.956471642313
E2F1	AACSL	1.0	3.6853160511
EGR1	AACSL	0.0	-0.695698519643

Example lioness output:
Sample1 Sample2 Sample3 Sample4
-------------------------------
-0.667452814003	-1.70433776179	-0.158129613892	-0.655795512803
-0.843366539284	-0.733709815256	-0.84849895139	-0.915217389738
3.23445386464	2.68888472802	3.35809757371	3.05297381396
2.39500370135	1.84608635425	2.80179804094	2.67540878165
-0.117475863987	0.494923925853	0.0518448588965	-0.0584810456421

TF, Gene and Motif order is identical to the panda output file.
```

# SAMBAR

# pySAMBAR #
## Subtyping Agglomerated Mutations By Annotation Relations ##

This package is the python implementation of the SAMBAR method as implemented in the R https://github.com/mararie/SAMBAR.
This package uses the python libraries numpy, pandas, networkx and scipy.

The easiest way to install pySAMBAR is by the following procedure using pip:
```
git clone https://github.com/genisott/pysambar.git
cd pysambar
pip install .
```
To use the package you can import it using: ```import pysambar```. And then access the different functions implemented with ```pysambar.function()```.

As an example you can find mutation data of Uterine Corpus Endometrial Carcinoma (UCEC) primary tumpor samples from The Cancer Genome Atlas. This data is in the ToyData folder as well as the MSigDb "Hallmark" gene sets. 

The program will compute the pathway mutation scores and clustering for *k*=2-4 (by default) and output the corrected mutation scores, the pathway mutation scores, and the clustering table. 

## More information on the method ##

SAMBAR, or **S**ubtyping **A**gglomerated **M**utations **B**y **A**nnotation **R**elations, is a method to identify subtypes based on somatic mutation data. SAMBAR was used to identify mutational subtypes in 23 cancer types from The Cancer Genome Atlas (Kuijjer ML, Paulson JN, Salzman P, Ding W, Quackenbush J, *British Journal of Cancer* (May 16, 2018), doi: 10.1038/s41416-018-0109-7, https://www.nature.com/articles/s41416-018-0109-7, *BioRxiv*, doi: https://doi.org/10.1101/228031).

SAMBAR's input is a matrix that includes the number of non-synonymous mutations in a sample <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;i" title="i" /></a> and gene <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;j" title="j" /></a>. SAMBAR first subsets these data to a set of 2,219 cancer-associated genes (optional) from the Catalogue Of Somatic Mutations In Cancer (COSMIC) and Östlund *et al*. (Network-based identification of novel cancer genes, 2010, *Mol Cell Prot*), or from a user-defined list. It then divides the number of non-synonymous mutations by the gene's length <a href="https://www.codecogs.com/eqnedit.php?latex=L_j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L_j" title="L_j" /></a>, defined as the number of non-overlapping exonic base pairs of a gene. For each sample, SAMBAR then calculates the overall cancer-associated mutation rate by summing mutation scores in all cancer-associated genes <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;j'" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;j'" title="j'" /></a>. It removes samples for which the mutation rate is zero and divides the mutation scores the remaining samples by the sample's mutation rate, resulting in a matrix of mutation rate-adjusted scores <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;G" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;G" title="G" /></a>:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;G_{ij}=\frac{N_{ij}/L_{j}}{\displaystyle\sum_{j'}({N_{ij'}/L_{j'}})}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;G_{ij}=\frac{N_{ij}/L_{j}}{\displaystyle\sum_{j'}({N_{ij'}/L_{j'}})}" title="G_{ij}=\frac{N_{ij}/L_{j}}{\displaystyle\sum_{j'}({N_{ij'}/L_{j'}})}" /></a>.

The next step in SAMBAR is de-sparsification of these gene mutation scores (agglomerated mutations) into pathway mutation (annotation relation) scores. SAMBAR converts a (user-defined) gene signature (.gmt format) into a binary matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;M" title="M" /></a>, with information of whether a gene <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;j" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;j" title="j" /></a> belongs to a pathway <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;q" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;q" title="q" /></a>. It then calculates pathway mutation scores <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;P" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;P" title="P" /></a> by correcting the sum of mutation scores of all genes in a pathway for the number of pathways <a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;q'" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;q'" title="q'" /></a> a gene belongs to, and for the number of cancer-associated genes present in that pathway:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{100}&space;P_{iq}=\frac{\displaystyle\sum_{j&space;\in&space;q}&space;G_{ij}/{\displaystyle\sum_{q'}&space;M_{jq'}}}{\displaystyle\sum_{j}&space;M_{jq}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dpi{100}&space;P_{iq}=\frac{\displaystyle\sum_{j&space;\in&space;q}&space;G_{ij}/{\displaystyle\sum_{q'}&space;M_{jq'}}}{\displaystyle\sum_{j}&space;M_{jq}}" title="P_{iq}=\frac{\displaystyle\sum_{j \in q} G_{ij}/{\displaystyle\sum_{q'} M_{jq'}}}{\displaystyle\sum_{j} M_{jq}}" /></a>.

Finally, SAMBAR uses binomial distance to cluster the pathway mutation scores. The cluster dendrogram is then divided into *k* groups (or a range of *k* groups), and the cluster assignments are returned in a list.
## Usage instructions ##
This package includes the functions ```sambar```,```desparsify```,```corgenelength```,```convertgmt```,```clustering``` as well as an implementation of the binomial distance (Millar dissimilarity from the package vegdist from R. To see the full description of each of this functions use ```help(pysambar.function)```.

## Example ## 

Run the SAMBAR method with the ToyData from the UCEC mutation data:
``` 
import pysambar as sm
pathways, groups = sm.sambar("/ToyData/mut.ucec.csv","ToyData/esizef.csv",'ToyData/genes.txt','ToyData/h.all.v6.1.symbols.gmt')
``` 
The output of this command will be four files:

``` pt_out.csv```  -> Pathway mutation score matrix.

``` mt_out.csv```  -> Processed gene mutation score matrix.

``` dist_matrix.csv ``` -> Distance matrix with binomial distance in numpy condensed format.

``` clustergroups.csv```  -> Matrix of pertinence to a cluster.

The function also returns the pathway matrix dataframe and the cluster group dataframe as python variables.
### Notes ###
Flags by default in the sambar function:

``` normPatient=True```  -> Normalizes the mutation data by number of mutations in a sample.

``` kmin=2,kmax=4```  -> Cut-offs of the cluster tree.

``` gmtMSigDB=True```  -> If the signature file comes from MSigDB the second element of each line is removed to process the format available at the database. This can be toggled off if using a custom signature file.

``` subcangenes=True```  -> Subset to cancer associated genes. By default uses the file provided in ToyData.

# CONDOR

# pyCONDOR
Python implementation of the BRIM algorithm for bipartite community structure detection as described in "Modularity and community detection in bipartite networks" by Michael J. Barber." This package is somewhat translated from the R package CONDOR https://github.com/jplatig/condor described in the paper "Bipartite Community Structure of eQTLs" by John Platig , Peter J. Castaldi, Dawn DeMeo, John Quackenbush.

## Install instructions
This package uses the python libraries numpy,pandas and python-igraph (igraph).

The easiest way to install pyCONDOR is using pip:
```
git clone https://github.com/genisott/pycondor.git
cd pycondor
pip install .
```
To use the package you can import it using: ```import condor```. And then access the different functions implemented with ```condor.function()```.


## Usage instructions
Supose you have a network (weighted or not) as an edgelist loaded into a pandas dataframe.
```
import condor
co = condor.condor_object(net)
```
Returns a condor object (dictionary with a graph, node names, community membership dataframes...)
```
co = condor.initial_community(co)
```
Computes the initial community structure and updates the condor object.
```
co = condor.brim(co)
```
Runs the iterative modularity optimization algorithm (BRIM) and updates the condor object with the final membership.
To see the results type:
```
co["tar_memb"]
co["reg_memb"]
```
To compute the qscores for the vertices type:
```
co = condor.qscores(co) # Computes the qscores
co["qscores"]["reg_qscores"] # Dataframe containing the qscores for the regulators.
co["qscores"]["tar_qscores"] # Dataframe containing the qscores for the targets.
```








