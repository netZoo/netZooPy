#!/usr/bin/env python

import sys
import getopt
import click
from netZooPy.panda.panda import Panda
from netZooPy.lioness import Lioness

#############################################################################
# PANDA #####################################################################
#############################################################################

@click.command()
@click.option('-e', '--expression', 'expression', type=str, required=True,
              help='Path to file containing the gene expression data. By default, \
                  the expression file does not have a header, and the cells are separated by a tab.')
@click.option('-m', '--motif', 'motif', type=str, required=True,
              help='Path to pair file containing the transcription factor DNA binding motif edges in the form of TF-gene-weight(0/1). If not provided, the gene coexpression matrix is returned as a result network.')
@click.option('-p', '--ppi', 'ppi', type=str, required=True,
              help='Path to pair file containing the PPI edges. The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.')
@click.option('-o', '--out', 'output', type=str, required=True,
              help='Output panda file. Format as txt')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--save_memory', is_flag=True, show_default=True,
              help='panda option. When true the result network is weighted adjacency matrix of size (nTFs, nGenes).\
                  when false The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.')
@click.option('--save_tmp', is_flag=True, show_default=True,
              help='panda option')
@click.option('--rm_missing', is_flag=True, show_default=False,
              help='Removes the genes and TFs that are not present in one of the priors. Works only if modeProcess=legacy')
@click.option('--keep_expr', is_flag=True, show_default=True,
              help='Keeps the input expression matrix in the result')
@click.option('--mode_process', type=str, default='union', show_default=True,
              help='panda option for input data processing. Choose between union(default), \
                  legacy and intersection')
@click.option('--alpha', type=float, default=0.1, show_default=True,
              help='panda and lioness first sample')
@click.option('--start', type=int, default=1, show_default=True,
              help='panda and lioness first sample')
@click.option('--end', type=int, default=None, show_default=True,
              help='panda and lioness last sample')
def panda(expression, motif, ppi, output, computing, precision, save_memory, save_tmp, rm_missing, keep_expr, mode_process, alpha, start, end):
    """Run panda using expression, motif and ppi data. 
    Use flags to modify the function behavior. By default, boolean flags are false.
    Output is a text file, with the TF, Gene, Motif, Force columns, where TF and Gene 
    are the nodes of the network considered, Motif is the prior and force is 
    the actual panda score.

    Example:

            netzoopy panda -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -o test_panda.txt
    
    Reference:

        Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PloS one 8.5 (2013): e64832.
   
    """
    print('Input data:')
    print('Expression:', expression)
    print('Motif data:', motif)
    print('PPI data:', ppi)

    # Run PANDA
    print('Start Panda run ...')
    panda_obj = Panda(expression, motif, ppi, computing=computing,precision = precision,  save_tmp=True, remove_missing=rm_missing, keep_expression_matrix=keep_expr, save_memory=save_memory,modeProcess=mode_process, alpha=alpha,start=start, end=end)
    panda_obj.save_panda_results(output)

#############################################################################
# LIONESS ###################################################################
#############################################################################

@click.command()
@click.option('-e', '--expression', 'expression', type=str, required=True,
              help='Path to file containing the gene expression data. By default, \
                  the expression file does not have a header, and the cells are separated by a tab.')
@click.option('-m', '--motif', 'motif', type=str, required=True,
              help='Path to pair file containing the transcription factor DNA binding motif edges in the form of TF-gene-weight(0/1). If not provided, the gene coexpression matrix is returned as a result network.')
@click.option('-p', '--ppi', 'ppi', type=str, required=True,
              help='Path to pair file containing the PPI edges. The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.')
@click.option('-op', '--out-panda', 'output_panda', type=str, required=True,
              help='Output panda file. Format as txt')
@click.option('-ol', '--out-lioness', 'output_lioness', type=str, required=True,
              help='Output lioness folder')
@click.option('--fmt', type=str, show_default=True, default='npy',
              help='Lioness network files output format. Choose one between .npy,.txt,.mat')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--ncores', type=int, show_default=True, default=1,
              help='Number of cores. Lioness CPU parallelizes over ncores')
@click.option('--save_memory', is_flag=True, show_default=True,
              help='panda option. When true the result network is weighted adjacency matrix of size (nTFs, nGenes).\
                  when false The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.')
@click.option('--save_tmp', is_flag=True, show_default=True,
              help='panda option')
@click.option('--rm_missing', is_flag=True, show_default=False,
              help='Removes the genes and TFs that are not present in one of the priors. Works only if modeProcess=legacy')
@click.option('--mode_process', type=str, default='union', show_default=True,
              help='panda option for input data processing. Choose between union(default), \
                  legacy and intersection')
@click.option('--output_type', type=str, default='network', show_default=True,
              help='lioness option for output format. Choose one between network, gene_targeting, tf_targeting')
@click.option('--alpha', type=float, default=0.1, show_default=True,
              help='panda and lioness first sample')
@click.option('--start', type=int, default=1, show_default=True,
              help='panda and lioness first sample')
@click.option('--end', type=int, default=None, show_default=True,
              help='panda and lioness last sample')
def lioness(expression, motif, ppi, output_panda, output_lioness, fmt, computing, precision, ncores, save_memory, save_tmp, rm_missing, mode_process,output_type, alpha, start, end):
    """Run Lioness to extract single-sample networks.
    First runs panda using expression, motif and ppi data. 
    Then runs lioness and puts results in the output_lioness folder.
    Use flags to modify the function behavior. By default, boolean flags are false.

    Example:

            netzoopy lioness -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -op test_panda.txt -ol lioness/
    
    Reference:
        Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.
    
    """
    
    print('Input data:')
    print('Expression:', expression)
    print('Motif data:', motif)
    print('PPI data:', ppi)
    print('Output panda:', output_panda)
    print('Lioness folder:', output_lioness)

    # Run PANDA
    print('Start Panda run ...')
    panda_obj = Panda(expression, motif, ppi, computing=computing, save_tmp=True, remove_missing=rm_missing, keep_expression_matrix=True, save_memory=save_memory, modeProcess=mode_process, start=start, end=end)
    panda_obj.save_panda_results(output_panda)
    print('Panda saved. Computing Lioness...')
    Lioness(panda_obj, computing=computing, precision=precision,ncores=ncores, save_dir=output_lioness, save_fmt=fmt, output = output_type, alpha = alpha)
    print('All done!')
    
