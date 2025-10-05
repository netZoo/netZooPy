#!/usr/bin/env python

import sys
import getopt
import click
from netZooPy.panda.panda import Panda
from netZooPy.lioness import Lioness
from netZooPy.bonobo import Bonobo
from netZooPy.condor import condor_object

#############################################################################
# PANDA #####################################################################
#############################################################################

def get_list_from_str(a):
    print(a)
    if a=='':
        return(None)
    else:
        return([int(b) for b in a.split(',')])

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
@click.option('--with_header', is_flag=True, show_default=False,
              help='Pass if the expression file has a header. It will be used to save samples with the correct name.')
@click.option('--save_memory', is_flag=True, show_default=True,
              help='panda option. When true the result network is weighted adjacency matrix of size (nTFs, nGenes).\
                  when false The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.')
@click.option('--as_adjacency', is_flag=True, show_default=True,
              help='If true, the final PANDA is saved as an adjacency matrix. Works only when save_memory is false')
@click.option('--old_compatible', is_flag=True, show_default=True,
              help='If true, PANDA is saved without headers. Pass this if you want the same results of netzoopy before v0.9.11')
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
              help='panda first sample')
@click.option('--end', type=int, default=None, show_default=True,
              help='panda last sample')
def panda(expression, motif, ppi, output, computing='cpu', precision='double',with_header=False, save_memory=False, as_adjacency=False, old_compatible=False, save_tmp=False, rm_missing=False, keep_expr=False, mode_process='union', alpha=0.1, start=1, end=None):
    """ Run panda using expression, motif and ppi data. 
    Use flags to modify the function behavior. By default, boolean flags are false.
    Output is a text file, with the TF, Gene, Motif, Force columns, where TF and Gene 
    are the nodes of the network considered, Motif is the prior and force is 
    the actual panda score.
    
    warning: To keep the command line call clean, we have set all booleans to false as default. To set a boolean to True,
    from command line you just have to add the respective flag
        >>> netzoopy panda ....  --rm_missing --save_tmp
    In this case rm_missing and save_tmp are set to True.
    To replicate the Panda class default behavior pass --save_memory and --save_tmp


    Example:

            netzoopy panda -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -o test_panda.txt


    Reference:

        Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PloS one 8.5 (2013): e64832.


    """
    print('NEW: We changed the default behavior of save_tmp (now False). Pass --save_tmp for retrocompatibility')
    print('NEW: We changed the default beavior of save_panda_results. Now all PANDA outputs have column headers by default. Pass old_compatible for the previous behavior')
    print('Input data:')
    print('Expression:', expression)
    print('Motif data:', motif)
    print('PPI data:', ppi)

    # Run PANDA
    print('Start Panda run ...')
    panda_obj = Panda(expression, motif, ppi, computing=computing,precision = precision,  save_tmp=save_tmp, remove_missing=rm_missing, keep_expression_matrix=keep_expr, save_memory=save_memory,modeProcess=mode_process, alpha=alpha,start=start, end=end, with_header=with_header)
    panda_obj.save_panda_results(output, save_adjacency=as_adjacency, old_compatible=old_compatible)

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
@click.option('--el', type=str, show_default=True, default='None',
              help='Lioness output export. If a file is passed, the final output will be saved as a complete table, with indices and column names, using the format specified here. Otherwise a standard liones.fmt file with no annotation is saved')
@click.option('--fmt', type=str, show_default=True, default='npy',
              help='Lioness network files output format. Choose one between .npy,.txt,.mat')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--ncores', type=int, show_default=True, default=1,
              help='Number of cores. Lioness CPU parallelizes over ncores')
#@click.option('--save_memory', is_flag=True, show_default=False,
#              help='panda option. When true the result network is weighted adjacency matrix of size (nTFs, nGenes).\
#                  when false The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.')
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
@click.option('--panda_start', type=int, default=1, show_default=True,
              help='panda first sample')
@click.option('--panda_end', type=int, default=None, show_default=True,
              help='panda last sample')
@click.option('--start', type=int, default=1, show_default=True,
              help='lioness first sample')
@click.option('--end', type=int, default=None, show_default=True,
              help='lioness last sample')
@click.option('--subset_numbers', type=str, default='', show_default=True,
              help='Specify a list of samples (numbers,comma separated) to run lioness on. \
                  Background is the one specified by panda_start and panda_end')
@click.option('--subset_names', type=str, default='', show_default=True,
              help='Specify a list of samples (sample names,comma separated) to run lioness on. \
                  Background is the one specified by panda_start and panda_end')
@click.option('--with_header', is_flag=True, show_default=False,
              help='Pass if the expression file has a header. It will be used to save samples with the correct name.')
@click.option('--save_single_lioness', is_flag=True, show_default=False,
              help='Pass this flag to save all single lioness networks generated.')
@click.option('--ignore_final', is_flag=True, show_default=False,
              help='The whole lioness data is not kept in memory. Always use save_single_lioness for this')
@click.option('--as_adjacency', is_flag=True, show_default=True,
              help='If true, the final PANDA is saved as an adjacency matrix. Works only when save_memory is false')
@click.option('--old_compatible', is_flag=True, show_default=True,
              help='If true, PANDA is saved without headers. Pass this if you want the same results of netzoopy before v0.9.11')

def lioness(expression, motif, ppi, output_panda, output_lioness, el, fmt, computing, precision, ncores, save_tmp, rm_missing, mode_process,output_type, alpha, panda_start, panda_end, start, end, subset_numbers='', subset_names='',with_header=False, save_single_lioness=False,ignore_final=False, as_adjacency=False, old_compatible=False):    

    """Run Lioness to extract single-sample networks.
    First runs panda using expression, motif and ppi data. 
    Then runs lioness and puts results in the output_lioness folder.
    Use flags to modify the function behavior. By default, boolean flags are false.

    Example:

            netzoopy lioness -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -op test_panda.txt -ol lioness/
    
    **Example for GPU computing**. Here the biggest limitation is GPU size, hence we need to optimize the precision and make
    sure we don't use all genes, but only the intersection with the PANDA priors. We also save the single lioness
    networks as they are computed:
    
        netzoopy lioness -e <expression_file> -m <motif_file> -p <ppi_file> -op output_panda.txt -ol output_lioness_folder/ --computing gpu --precision single --mode_process intersection --save_single_lioness --ignore_final

    **LIONESS on a subset of samples**. This is especially useful if you need to test whether your data is suitable for
    LIONESS/you have enough resources. By specifying --panda_start and --panda_end the number of samples are restricted
    and only the samples in the subset are used for LIONESS. The background is the one specified by panda_start and
    panda_end. Alternatively you can use the --subset_numbers and --subset_names flags to specify the samples to use. In this case
    the PANDA will be computed on all samples, but then LIONESS will be computed and saved only for a subset of samples:
    
        netzoopy lioness -e <expression_file> -m <motif_file> -p <ppi_file> -op output_panda.txt -ol output_lioness_folder/ --computing <gpu|cpu> --panda_start 1 --panda_end 5 --precision single --mode_process intersection --save_single_lioness --ignore_final
    


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
    
    # For now we keep save_memory=False always, otherwise we won't have the needed information for lioness
    panda_obj = Panda(expression, motif, ppi, precision=precision, computing=computing, save_tmp=save_tmp, remove_missing=rm_missing, keep_expression_matrix=True, save_memory=False, modeProcess=mode_process, start=panda_start, end=panda_end, with_header=with_header)
    print('Panda saved. Computing Lioness...')
    panda_obj.save_panda_results(output_panda, save_adjacency=as_adjacency, old_compatible=old_compatible)

    if el=='None':
        el = None

    subset_numbers = get_list_from_str(subset_numbers)
    subset_names = get_list_from_str(subset_names)

    Lioness(panda_obj, computing=computing, precision=precision,ncores=ncores, save_dir=output_lioness, save_fmt=fmt, output = output_type, alpha = alpha, export_filename=el, save_single=save_single_lioness,ignore_final=ignore_final, start=start, end=end, subset_names=subset_names, subset_numbers=subset_numbers)

    print('All done!')
    

################################################
###### CONDOR ##################################
################################################

@click.command()
@click.option('-n', '--network_file', 'network', type=str, required=True,
              help='Path to file encoding an edgelist.')
@click.option('--sep', type=str, show_default=True, default=',',
              help='network file separator')
@click.option('--index_col', type=int, show_default=True, default=0,
              help='Column that stores the index of the edgelist. E.g. None, 0...')
@click.option('--header', type=int, show_default=True, default=0,
              help='Row that stores the header of the edgelist. E.g. None, 0...')
@click.option('--initial_method', type=str, show_default=True, default='LCS',
              help='Method to determine intial community assignment. (By default Multilevel method)')
@click.option('--initial_project', is_flag=True, show_default=True,
              help='Whether to project the network onto one of the bipartite sets for\
                    the initial community detection.')       
@click.option('--com_num', type=str, show_default=True, default='def',
              help='Max number of communities. It is recomended to leave this to default,\
                   otherwise if the initial community assignement is bigger the program will crash.')
@click.option('--delta_qmin', type=str, show_default=True, default='def',
              help='Difference modularity threshold for stopping the iterative process.')
@click.option('--resolution', type=int, show_default=True, default=1,
              help='Not yet implemented')
@click.option('--tar_output', type=str, show_default=True, default='tar_memb.txt',
              help='Filename for saving the tar node final membership.')
@click.option('--reg_output', type=str, show_default=True, default='reg_memb.txt',
              help='Filename for saving the reg node final membership.')
def condor(
    network_file,
    sep=",",
    index_col=0,
    header=0,
    initial_method="LCS",
    initial_project=False,
    com_num="def",
    delta_qmin="def",
    resolution=1,
    tar_output="tar_memb.txt",
    reg_output="reg_memb.txt",
):
    """
        Computation of the whole condor process. It creates a condor object and runs all the steps of BRIM on it. The function outputs
        Note: The edgelist is assumed to contain a bipartite network. The program will relabel the nodes so that the edgelist represents a bipartite network anyway.
        It is on the user to know that the network they are using is suitable for the method.
        
        tar and reg final memberships are saved to csv
    """

    co = condor_object(network_file, sep, index_col, header)
    co.initial_community(initial_method, initial_project)
    co.brim(delta_qmin, com_num, resolution)
    co.tar_memb.to_csv(tar_output)
    co.reg_memb.to_csv(reg_output)



################################################
###### BONOBO ##################################
################################################

@click.command()
@click.option('-e', '--expression_file', 'expression_file', type=str, required=True,
              help='Path to file containing the gene expression data or pandas dataframe. By default, the expression file does not have a header, and the cells ares separated by a tab.')
@click.option('--output_folder', type=str, show_default=True, default='bonobo/',
              help='Output folder for the bonobo files. If not specified, the bonobo files are saved in the current directory, in the bonobo subdirectory.')
@click.option('--output_format', type=str, show_default=True, default='.h5',
              help='format of output bonobo matrix. By default it is an hdf file, can be a txt or csv.')
@click.option('--keep_in_memory', is_flag=True, show_default=True,
              help='if True, the bonobo coexpression matrix is kept in memory, otherwise it is discarded after saving')    
@click.option('--delta', type=float, show_default=True, default=None,
              help='delta parameter. If default (None) delta is trained, otherwise pass a value.Recommended is 0.3.')
@click.option('--sparsify', is_flag=True, show_default=True,
              help='if True, bonobo gets sparsified and relative pvalues are returned')    
@click.option('--confidence', type=float, show_default=True, default=0.05,
              help='if sparsify is True, this is the CI for the approximate zscore.')
@click.option('--save_pvals', is_flag=True, show_default=True,
              help='if True, bonobo gets sparsified and relative pvalues are saved in the same format and folder of bonobo')    
@click.option('--precision', type=str, show_default=True, default='single',
              help='matrix precision (single or double), defaults to single precision.')
@click.option('--sample_names', type=str, show_default=True, default='',
              help='Compute BONOBO only on a subset of samples. Pass a comma separated list of sample names. If not specified, all samples are used.')
def bonobo(
    expression_file,
    output_folder = 'bonobo/',
    output_format = '.h5',
    keep_in_memory = False, 
    delta = None, 
    sparsify = False, 
    confidence = 0.05, 
    save_pvals = False,
    precision = 'single',
    sample_names = '',
):
    """
        Compute BONOBOs from an expression file. 
        
        Parameters the user cannot access from the CLI:
        - computing: for now it is only CPU
        - cores: number of cores to use, for now there is no parallelization
        - online_coexpression: we have not implemented the online coexpression yet
    """

    if sample_names!='':
        sample_names = sample_names.split(',')
        print('WARNING: computing BOBOBO only on a subset of samples. The sample names are:')
        print(sample_names)
    else:
        sample_names = []

    print('Initializing BONOBO object ...')
    bonobo_obj_sparse = Bonobo(expression_file)
    print('Running BONOBO ...')
    print('Files saved in %s' %output_folder)

    bonobo_obj_sparse.run_bonobo(keep_in_memory=keep_in_memory, 
                                 output_fmt=output_format, 
                                 delta = delta, 
                                 sparsify=sparsify, 
                                 output_folder=output_folder, 
                                 confidence = confidence,
                                 save_pvals=save_pvals, 
                                 precision = precision, 
                                 sample_names=sample_names)



#####################################################################################
############## OTTER LIONESS ########################################################
#####################################################################################
    
from netZooPy.lioness.lioness_for_otter import LionessOtter

@click.command()
@click.option('-e', '--expression', 'expression', type=str, required=True,
              help='Path to file containing the gene expression data. By default, \
                  the expression file does not have a header, and the cells are separated by a tab.')
@click.option('-m', '--motif', 'motif', type=str, required=True,
              help='Path to pair file containing the transcription factor DNA binding motif edges in the form of TF-gene-weight(0/1). If not provided, the gene coexpression matrix is returned as a result network.')
@click.option('-p', '--ppi', 'ppi', type=str, required=True,
              help='Path to pair file containing the PPI edges. The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.')
@click.option('-of', '--out-folder', 'output_folder', type=str, required=True,
              help='Output lioness otter folder')
@click.option('--fmt', type=str, show_default=True, default='h5',
              help='Lioness network files output format. Choose one between .npy,.txt,.mat')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--mode_process', type=str, default='intersection', show_default=True,
              help='panda option for input data processing. Choose between union(default), \
                  legacy and intersection')
@click.option('--iterations', type=int, default=60, show_default=True,
              help='otter iterations, Iter')
@click.option('--lam', type=float, default=0.035, show_default=True,
              help='lambda parameter')
@click.option('--gamma', type=float, default=0.335, show_default=True,
              help='gamma parameter')
@click.option('--eta', type=float, default=0.00001, show_default=True,
              help='eta parameter')
@click.option('--bexp', type=int, default=1., show_default=True,
              help='bexp parameter')
def otterlioness(expression, motif, ppi, output_folder, fmt, computing, precision, mode_process='intersection', iterations=60, lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
    """Run Lioness otter to extract single-sample networks.
    First runs otter using expression, motif and ppi data. 
    Then runs lioness and puts results in the output_lioness folder.
    WARNING: the OTTER CLI and class are still relying on a simple approach for reading and merging. Please be careful
    if you have NAs and want a non-intersection between W,P,C please rely on PANDA or on your own filtering. 

    Example:

            netzoopy otterlioness -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -of lioness_otter/
    
    """
    # Run PANDA
    print('Start Otter run ...')

    # First we create the LIONESS OTTER instance with the expression, motif, ppi files        
    lioobj = LionessOtter(expression, motif, ppi, mode_process=mode_process)
    
    print('Starting Otter Lioness')
    lioobj.run_lioness_otter(output_folder, save_fmt = fmt, save_single=True, precision = precision, computing = computing, Iter = iterations, lam=lam, gamma=gamma, eta=eta, bexp=bexp)




#####################################################################################
############## OTTER ################################################################
#####################################################################################
    
from netZooPy.lioness.lioness_for_otter import LionessOtter

@click.command()
@click.option('-e', '--expression', 'expression', type=str, required=True,
              help='Path to file containing the gene expression data. By default, \
                  the expression file does not have a header, and the cells are separated by a tab.')
@click.option('-m', '--motif', 'motif', type=str, required=True,
              help='Path to pair file containing the transcription factor DNA binding motif edges in the form of TF-gene-weight(0/1). If not provided, the gene coexpression matrix is returned as a result network.')
@click.option('-p', '--ppi', 'ppi', type=str, required=True,
              help='Path to pair file containing the PPI edges. The PPI can be symmetrical, if not, it will be transformed into a symmetrical adjacency matrix.')
@click.option('-o', '--out-file', 'output_file', type=str, default = 'otter.txt',
              help='Output otter file. Use one of the extensions between .npy,.txt,.mat')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--mode_process', type=str, default='intersection', show_default=True,
              help='panda option for input data processing. Choose between union(default), \
                  legacy and intersection')
@click.option('--iterations', type=int, default=60, show_default=True,
              help='otter iterations, Iter')
@click.option('--lam', type=float, default=0.035, show_default=True,
              help='lambda parameter')
@click.option('--gamma', type=float, default=0.335, show_default=True,
              help='gamma parameter')
@click.option('--eta', type=float, default=0.00001, show_default=True,
              help='eta parameter')
@click.option('--bexp', type=int, default=1., show_default=True,
              help='bexp parameter')
def otter(expression, motif, ppi, output_file='otter.txt', computing='cpu', precision='double', mode_process='intersection', iterations=60, lam=0.035, gamma=0.335, Iter=60, eta=0.00001, bexp=1):
    """Run Lioness otter to extract single-sample networks.
    First runs otter using expression, motif and ppi data. 
    Then runs lioness and puts results in the output_lioness folder.
    
    WARNING: the OTTER CLI and class are still relying on a simple approach for reading and merging. Please be careful
    if you have NAs and want a non-intersection between W,P,C please rely on PANDA or on your own filtering. 
    Example:

            netzoopy otterlioness -e tests/puma/ToyData/ToyExpressionData.txt -m tests/puma/ToyData/ToyMotifData.txt -p tests/puma/ToyData/ToyPPIData.txt -of lioness_otter/
    
    """
    # Run PANDA
    print('Start Otter run ...')

    # First we create the LIONESS OTTER instance with the expression, motif, ppi files        
    lioobj = LionessOtter(expression, motif, ppi, mode_process=mode_process)
    
    print('Starting Otter Lioness')
    
    lioobj.run_otter(output_file, precision = precision, computing = computing, Iter = iterations, lam=lam, gamma=gamma, eta=eta, bexp=bexp )

