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










