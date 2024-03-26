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

def check_expression_integrity(df):
    """Check data integrity
    - Number of NA

    Args:
        df (dataframe): gene expression dataframe
    """

    # check that for each
    if (df.isna().sum(axis = 1)>(len(df.columns)-3)).any():
        sys.exit('Too many nan in gene expression (need more than 1 sample to compute coexpression)')

def read_ppi(ppi_fn, tf_list = None):
    """Read PPI network

    Args:
        ppi_fn (str): ppi network filename
    """
    with open(ppi_fn, 'r') as f:
        ppi_data = pd.read_csv(f, sep="\t", header=None)
        ppi_data.columns = ['tf1','tf2','exists']

    # get all tfs from first and second column
    if tf_list:
        ppi_tfs = tf_list
        ppi_data = ppi_data[(ppi_data.tf1.isin(ppi_tfs)) & (ppi_data.tf2.isin(ppi_tfs))]
    else:
        ppi_tfs = sorted(set(ppi_data.iloc[:,0].values.tolist()).union(set(ppi_data.iloc[:,1].values.tolist())))
    
    # create adjacency matrix
    df = pd.DataFrame(np.eye(len(ppi_tfs)), index=ppi_tfs, columns=ppi_tfs)
    z = ppi_data.pivot_table(columns='tf2',index = 'tf1',values = 'exists', fill_value=0)
    df = df.add(z, fill_value=0).add(z.T, fill_value=0)
    df = 1*(df>0)
    # return adjacency matrix and tfs list
    return(df, ppi_tfs)

def read_motif(motif_fn, pivot = True):
    """ Read a motif edgelist, generates

    Args:
        motif_fn (_type_): filename of the motif edgelist
        pivot (bool): if true returns a pivot tfs X genes table. Otherwise keeps the edgelist
    Returns:
        piv/df: motif as edgelist or pivot table
        tfs: list of tfs
        genes: list of genes
    """

    with open(motif_fn, 'r') as f:
        df = pd.read_csv(f, sep= '\t', header = None)

    presenttf = df.iloc[:,0].unique()
    presentgene = df.iloc[:,1].unique()

    if pivot:
        piv = df.pivot_table(values=2, index=0, columns=1, fill_value=0)
        return(piv, list(presenttf), list(presentgene))
    else:
        return(df, list(presenttf), list(presentgene))

def read_expression(expression_fn, header = 0, usecols = None, nrows = None):
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
            df = pd.read_csv(f, sep = '\t', usecols = usecols, index_col=0, nrows=nrows)
        elif expression_fn.endswith('.csv'):
            df = pd.read_csv(f, sep = ' ', usecols = usecols, index_col=0, nrows=nrows)
        elif expression_fn.endswith('.tsv'):
            df = pd.read_csv(f, sep = '\t', usecols = usecols, index_col=0, nrows=nrows)
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
    print(samples)
    if type(expression_filename) is str:
        columns = read_expression(expression_filename, nrows = 1)
        if (isinstance(samples, list)):
            usecols = samples.copy()
            usecols.insert(0,columns.index.name)
            expression_data = read_expression(expression_filename, usecols = usecols)
        else:
            expression_data = read_expression(expression_filename)

    elif isinstance(expression_filename, pd.DataFrame):
        if (isinstance(samples, list)):
            usecols = samples.copy()
            usecols.insert(0,columns.index.name)
            expression_data = expression_filename.loc[:,usecols]
        else:
            expression_data = expression_filename
    else: 
        sys.exit('Expression filename needs to be either a table string or a panda dataframe')
    
    # keep names of expression genes
    expression_genes = set(expression_data.index.tolist())

    if len(expression_data) != len(expression_genes):
        print(
            "Duplicate gene symbols detected. Consider averaging before running PANDA"
        )

    check_expression_integrity(expression_data)

    return(expression_data, expression_genes)