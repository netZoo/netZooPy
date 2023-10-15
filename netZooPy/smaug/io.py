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


def read_expression(expression_fn, header=0, usecols=None, nrows=None):
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
            df = pd.read_csv(f, sep='\t', usecols=usecols, index_col=0, nrows=nrows)
        elif expression_fn.endswith('.csv'):
            df = pd.read_csv(f, sep=' ', usecols=usecols, index_col=0, nrows=nrows)
        elif expression_fn.endswith('.tsv'):
            df = pd.read_csv(f, sep='\t', usecols=usecols, index_col=0, nrows=nrows)
        else:
            sys.exit("Format of expression filename not recognised %s" % str(expression_fn))

    return (df)

def prepare_expression(expression_filename, samples=None):
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

        if (isinstance(samples, list)):
            expression_data = read_expression(expression_filename, usecols=samples)
        else:
            expression_data = read_expression(expression_filename)

    elif isinstance(expression_filename, pd.DataFrame):
        if (isinstance(samples, list)):
            expression_data = expression_filename.loc[:, samples]
        else:
            expression_data = expression_filename

    else:
        sys.exit('Filename needs to be a table string')

    # keep names of expression genes
    expression_genes = set(expression_data.index.tolist())

    if len(expression_data) != len(expression_genes):
        print(
            "Duplicate symbols detected. Consider averaging before running SMAUG"
        )

    check_expression_integrity(expression_data)

    return (expression_data, expression_genes)
