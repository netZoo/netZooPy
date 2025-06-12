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

def check_data_integrity(df):
    """Check data integrity
    - Number of NA

    Args:
        df (dataframe): expression or methylation dataframe
    """

    # check that for each
    if (df.isna().sum(axis = 1)>(len(df.columns)-3)).any():
        sys.exit('Too many nan in data (need more than 1 sample to compute partial correlation)')


def read_data(data_file_name, header=0, usecols=None, nrows=None):
    """Read data.

    Parameters
    -----------
        data_file_name: str
            filename of the expression or methylation file
        header: str or int
            header row
        usecols:list
            pass a list of the columns that need to be read
    """
    with open(data_file_name, 'r') as f:
        if data_file_name.endswith('.txt'):
            df = pd.read_csv(f, sep='\t', usecols=usecols, index_col=0, nrows=nrows)
        elif data_file_name.endswith('.csv'):
            df = pd.read_csv(f, sep=' ', usecols=usecols, index_col=0, nrows=nrows)
        elif data_file_name.endswith('.tsv'):
            df = pd.read_csv(f, sep='\t', usecols=usecols, index_col=0, nrows=nrows)
        else:
            sys.exit("Format of data filename not recognised %s" % str(expression_fn))

    return (df)

def prepare_data(data_filename, samples=None):
    """ Prepare main coexpression network by reading the expression file.

    Parameters
    ----------
        data_filename :str
            A table (tsv, csv, or txt) where each column is a sample
            and each row is a gene. Values are expression or methylation M values.
        samples: list
            list of sample names. If None all samples are read (default: None)

	Returns
	---------
		omics_data: pd.DataFrame
		omics_genes:set

    """
    # expression file is properly annotated with the sample name and
    # a list of sample of interest is passed
    print(samples)
    if type(data_filename) is str:

        if (isinstance(samples, list)):
            omics_data = read_data(data_filename, usecols=samples)
        else:
            omics_data = read_data(data_filename)

    elif isinstance(data_filename, pd.DataFrame):
        if (isinstance(samples, list)):
            omics_data = data_filename.loc[:, samples]
        else:
            omics_data = data_filename

    else:
        sys.exit('Filename needs to be a table or string')

    # keep names of expression genes
    omics_genes = set(omics_data.index.tolist())

    if len(omics_data) != len(omics_genes):
        print(
            "Duplicate symbols detected. Consider averaging before running SMAUG"
        )

    check_data_integrity(omics_data)

    return (omics_data, omics_genes)
