from __future__ import print_function
import math
from random import sample
import time
import pandas as pd
from scipy.stats import zscore
from .timer import Timer
import numpy as np
from netZooPy.panda import calculations as calc
import sys


def read_ppi(ppi_fn):

    with open(ppi_fn, 'r') as f:
        ppi_data = pd.read_csv(f, sep="\t", header=None)
        ppi_data.columns = ['tf1','tf2','exists']

    # get all tfs from first and second column
    ppi_tfs = sorted(set(ppi_data.iloc[:,0].values.tolist()).union(set(ppi_data.iloc[:,1].values.tolist())))
    # create adjacency matrix
    df = pd.DataFrame(np.eye(len(ppi_tfs)), index=ppi_tfs, columns=ppi_tfs)
    z = ppi_data.pivot_table(columns='tf2',index = 'tf1',values = 'exists', fill_value=0)
    df = df.add(z, fill_value=0).add(z.T, fill_value=0)
    df = 1*(df>0)
    return(df, ppi_tfs)


def read_priors_table(table_filename, sample_col = 'sample', prior_col = 'prior'):
    """Read the table detailing the sample-prior dyads

    Parameters
    ----------
        table_filename :str
            A csv that couples each sample with a prior
            filename. The csv file can have any number of columns, but only 
            the sample_col and prior_col are read. 
        sample_col: str 
            name of the column storing the sample name (default: 'sample')
        prior_col: str
            name of the column storing the prior filename  (default: 'prior')

    Returns
    -------
        samples: list
            list of samples that are defined
        prior_dict: dict
            dictionary sample:prior_filename
    """    

    with open(table_filename, 'r') as f:
        df = pd.read_csv(table_filename, usecols=[sample_col, prior_col])
    
    # check that sample-prior are unique
    if len(df.drop_duplicates(subset = [sample_col, prior_col]))!=len(df):
        print('Some samples are reported twice')
        print(df[df.duplicated(subset = [sample_col, prior_col])])
        df = df.drop_duplicates(subset = [sample_col, prior_col])

    # check that there are no more samples duplicated
    if len(df.drop_duplicates(subset = [sample_col]))!=len(df):
        sys.exit('No unique sample-prior assignment, there are duplicated sample columns')

    samples = df[sample_col].astype(str).values.tolist()
    sample2prior_dict = {i[1]:i[2] for i in df.itertuples()}

    prior2sample_dict = {}
    for p,tab in df.groupby(prior_col):
        prior2sample_dict[p] = tab.loc[:,sample_col].values.tolist()

    return(samples,  sample2prior_dict,  prior2sample_dict)

def read_motif(motif_fn, tf_names = None, gene_names = None, pivot = True):
    """ Read a motif edgelist, generates

    Args:
        motif_fn (_type_): filename of the motif edgelist
        tf_names (_type_, optional): list of tf_names to be used. Defaults to None.
        gene_names (_type_, optional): list of gene_names to be used. Defaults to None.
        pivot (bool): if true returns a pivot tfs X genes table. Otherwise keeps the edgelist
    """

    with open(motif_fn, 'r') as f:
        df = pd.read_csv(f, sep= '\t', header = None)

    tftoadd = None
    genetoadd = None

    presenttf = df.iloc[:,0].unique()
    presentgene = df.iloc[:,1].unique()
    ntfs = len(presenttf)
    ngenes = len(presentgene)

    if isinstance(tf_names, list):
        
        # Check how many are removed
        dtf = len(set(presenttf).difference(set(tf_names)))
        if (dtf!=0):
            print('Note: %d/%d tfs in the prior are not in the tf universe ' %(int(dtf),int(ntfs)))

        # tfs to add
        tftoadd = set(tf_names).difference(set(presenttf))

    if isinstance(gene_names, list):
    

        # Check how many are removed
        dgene = len(set(presentgene).difference(set(gene_names)))
        if (dgene!=0):
            print('Note: %d/%d genes in the prior are not in the gene universe ' %(int(dgene), int(ngenes)))

        # genes to add
        genetoadd = set(gene_names).difference(set(presentgene))

    print('hello')
    # now if one or both tftoadd/genetoadd are not None, we add rows with zeros at the end
    if (tftoadd or genetoadd):
        print('here')
        # if one of the two is None, add all the existing ones
        if tftoadd==None:
            tftoadd = presenttf
        if genetoadd==None:
            genetoadd = presentgene

        toadd = pd.DataFrame(tftoadd, columns = [0]).merge(pd.DataFrame(genetoadd, columns = [1]), how='cross')
        toadd[2] = 0

        df = pd.concat([df,toadd])

    if pivot:
        piv = df.pivot_table(values=2, index=0, columns=1, fill_value=0)
        if isinstance(tf_names, list):
            piv = piv.loc[tf_names,:]
        if isinstance(gene_names, list):
            piv = piv.loc[:,gene_names]
        return(piv)
    else:
        if isinstance(tf_names, list):
            df = df[df[0].isin(tf_names)]
        if isinstance(gene_names, list):
            df = df[df[1].isin(gene_names)]
        return(df)

def read_expression(expression_fn, header = 0, usecols = None):
    """Read expression data.

    Parameters
    -----------
        expression_fn: str
            filename of the expression file
        header: str or int
            header row
        usecols:list
            pass a list of the columns that need to be read
    """
    with open(expression_fn, 'r') as f:
        if expression_fn.endswith('.txt'):
            df = pd.read_csv(f, sep = '\t', usecols = usecols, index_col=0)
        elif expression_fn.endswith('.csv'):
            df = pd.read_csv(f, sep = ' ', usecols = usecols, index_col=0)
        elif expression_fn.endswith('.tsv'):
            df = pd.read_csv(f, sep = '\t', usecols = usecols, index_col=0)
        else:
            sys.exit("Format of expression filename not recognised %s" %str(expression_fn))
    
    return(df)


def prepare_expression(expression_filename, samples = None):

    """ Prepare main coexpression network by reading the expression file.
    
    Parameters
    ----------
        expression_filename :str
            A table (tsv, csv, or txt) where each column is a sample 
            and each row is a gene. Values are expression.
        samples: list
            list of sample names. If None all samples are read (default: None)

	Returns
	---------
		expression_data: pd.DataFrame
		expression_genes:set

    """    
    # expression file is properly annotated with the sample name and 
    # a list of sample of interest is passed
    
    if type(expression_filename) is str:
        if (isinstance(sample, list)):
            expression_data = read_expression(expression_filename, usecols = samples)
        else:
            expression_data = read_expression(expression_filename)

        
    elif isinstance(expression_filename, pd.DataFrame):

        expression_data = expression_filename.loc[:,samples]

    else: 
        sys.exit('Expression filename needs to be either a table string or a panda dataframe')
    

    # keep names of expression genes
    expression_genes = set(expression_data.index.tolist())

    if len(expression_data) != len(expression_genes):
        print(
            "Duplicate gene symbols detected. Consider averaging before running PANDA"
        )

    return(expression_data, expression_genes)

def read_motif_universe(priors_dict, mode = 'union', tf_col = 0, gene_col = 1):

    """ Read tf and gene names for all the priors and generate 
    a universe. By default the union of all gene names and tfs 
    is considered. Tables needs to be tab separated.

    Parameters
    ----------- 
    priors_dict: dict
        dictionary sample:prior_filename
    mode: str
        how to generate the universe. Default: union
	tf_col:str or int
		Column used to store the TF names
	gene_col: str or int
		Column used to store the gene names

    Returns
    --------
    tfs:set
        names of tfs
	genes:set
		names of genes
	universe:set
		name of both genes and tfs
    """

    files = set(priors_dict.values())
    tfs = set()
    genes = set()
    for fff in files:
        with open(fff, 'r') as f:
            df = pd.read_csv(f, sep= '\t', header = None)
            tf = df.iloc[:,0].unique().tolist()
            gene = df.iloc[:,1].unique().tolist()
            if mode=='union':
                tfs = tfs.union(set(tf))
                genes = genes.union(set(gene))
                universe = genes.union(tfs)
            elif mode=='intersection':
                tfs = tfs.intersection(set(tf))
                genes = genes.intersection(set(gene))   
                universe = genes.intersection(tfs)       
            else:
                sys.exit('Name %s is not a valid mode for prior tf/gene universe. Use union or intersection' %str(mode))   
    return(tfs, genes)




