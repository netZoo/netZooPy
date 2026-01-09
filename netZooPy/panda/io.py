from __future__ import print_function

import math
import numpy as np
import sys
import pandas as pd


def load_motif(motif_file):
    if type(motif_file) is str:
        # If motif_file is a filename
        motif_data = pd.read_csv(motif_file, sep="\t", header=None)
        motif_tfs = sorted(set(motif_data[0]))
        motif_genes = sorted(set(motif_data[1]))
    elif type(motif_file) is not str:
        # If motif_file is an object
        if motif_file is None:
            # Computation without motif
            motif_data = None
            motif_genes = []
            motif_tfs = []
        else:
            # If motif_file is an object, it needs to be a dataframe
            if not isinstance(motif_file, pd.DataFrame):
                raise Exception(
                    "Please provide a pandas dataframe for motif data with column names as 'source', 'target', and 'weight'."
                )
            if ("source" not in motif_file.columns) or (
                "target" not in motif_file.columns
            ):
                print('renaming motif columns to "source", "target" and "weight" ')
                motif_file.columns = ["source", "target", "weight"]
            motif_data = pd.DataFrame(motif_file.values)
            motif_tfs = sorted(set(motif_file["source"]))
            motif_genes = sorted(set(motif_file["target"]))
    return(motif_data, motif_tfs, motif_genes)


def load_expression(expression_file, with_header = False, start = 1, end = None):
    if type(expression_file) is str:
        # If we pass an expression file, check if we have a 'with header' flag and read it
        
            if with_header:
                # Read data with header
                expression_data = pd.read_csv(
                    expression_file, sep="\t", index_col=0
                )
            else:
                expression_data = pd.read_csv(
                    expression_file, sep="\t", header=None, index_col=0
                )
            # assign expression data and samples/gene names
            expression_data = expression_data.iloc[:, (start-1):end]
            expression_genes = expression_data.index.tolist()
            expression_samples = expression_data.columns.astype(str)
            
    elif type(expression_file) is not str:
        # Pass expression as a dataframe 
        if expression_file is not None:
            if not isinstance(expression_file, pd.DataFrame):
                raise Exception(
                    "Please provide a pandas dataframe for expression data."
                )
            expression_data = expression_file.iloc[:, (start-1):end]  # pd.read_csv(expression_file, sep='\t', header=None, index_col=0)
            expression_genes = expression_data.index.tolist()
            expression_samples = expression_data.columns.astype(str)
        else:
            expression_data = None
            expression_genes = None
            expression_samples = None
            print('No valid expression file is passed')

    return(expression_data, expression_genes, expression_samples)


#def load_ppi(ppi_file):
