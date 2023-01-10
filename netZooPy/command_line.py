#!/usr/bin/env python

import sys
import getopt
import click
from netZooPy.panda.panda import Panda
from netZooPy.lioness import Lioness
from netZooPy.ligress import Ligress
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
@click.option('--save_memory', is_flag=True, show_default=False,
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
def lioness(expression, motif, ppi, output_panda, output_lioness, el, fmt, computing, precision, ncores, save_memory, save_tmp, rm_missing, mode_process,output_type, alpha, panda_start, panda_end, start, end, subset_numbers='', subset_names='',with_header=False, save_single_lioness=False,ignore_final=False, as_adjacency=False, old_compatible=False):
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
    
    panda_obj = Panda(expression, motif, ppi, precision=precision, computing=computing, save_tmp=save_tmp, remove_missing=rm_missing, keep_expression_matrix=True, save_memory=save_memory, modeProcess=mode_process, start=panda_start, end=panda_end, with_header=with_header)
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



#############################################################################
# TIGRESS ###################################################################
#############################################################################

@click.command()
@click.option('-e', '--expression', 'expression', type=str, required=True,
              help='Path to file containing the gene expression data. By default, \
                  the expression file does not have a header, and the cells are separated by a tab.')
@click.option('-pt', '--priors-table', 'priors_table', type=str, required=True,
              help='Path to file containing the sample-prior table.' )
@click.option('-ol', '--out-lioness', 'output_lioness', type=str, required=True,
              help='Output lioness folder')
@click.option('-p','--ppi', type=str, show_default=True, default=None,
              help='PPI table file. Sample-PPI table')
@click.option('--ppitable', type=str, show_default=True, default='',
              help='PPI table file. Sample-PPI table')
@click.option('--fmt', type=str, show_default=True, default='npy',
              help='Lioness network files output format. Choose one between .npy,.txt,.mat')
@click.option('--computing', type=str, show_default=True, default='cpu',
              help='computing option, choose one between cpu and gpu')
@click.option('--precision', type=str, show_default=True, default='double',
              help='precision option')
@click.option('--ncores', type=int, show_default=True, default=1,
              help='Number of cores. Lioness CPU parallelizes over ncores')
@click.option('--save_memory', is_flag=True, show_default=False,
              help='panda option. When true the result network is weighted adjacency matrix of size (nTFs, nGenes).\
                  when false The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge.')
@click.option('--save_coexpression', is_flag=True, show_default=True,
              help='save ')
@click.option('--rm_missing', is_flag=True, show_default=False,
              help='Removes the genes and TFs that are not present in one of the priors. Works only if modeProcess=legacy')
@click.option('--mode_process', type=str, default='union', show_default=True,
              help='panda option for input data processing. Choose between union(default), \
                  legacy and intersection')
@click.option('--alpha', type=float, default=0.1, show_default=True,
              help='panda and lioness first sample')
@click.option('--output_type', type=str, default='txt', show_default=True,
              help='Output type. Now unused')
@click.option('--th_motifs', type=int, show_default=True, default=3,
              help='Threshold for the motifs. Reads motif only once if possible.')      
@click.option('--delta', type=float, default=0.3, show_default=True,
              help='posterior weight for single-sample coexpression estimation')
@click.option('--tune_delta', is_flag=True, show_default=True,
              help='tune the posterior weight for the single-sample coexpression estimation')
def ligress(expression, priors_table, ppi,output_lioness,ppitable, fmt, computing, precision, ncores, save_memory, save_coexpression, rm_missing, mode_process,output_type, alpha, th_motifs, delta, tune_delta):
    """Run Lioness to extract single-sample coexpression networks. 
    Then run Panda on each sample with sample-specific priors.
    """
    
    print('Input data:')
    print('Expression:', expression)

    # read expression data, prepare ppi+motif+expression universes
    print('Initialise ligress...')
    if ppitable!="":
        ligress_obj = Ligress(expression, priors_table, ppi_table_file=ppitable, output_folder=output_lioness, mode_process=mode_process)
    else:
        ligress_obj = Ligress(expression, priors_table, ppi_file = ppi, output_folder=output_lioness, mode_process=mode_process)
     
    print('Running ligress computations ...')
    ligress_obj.run_ligress(keep_coexpression=save_coexpression,cores=ncores, save_memory=save_memory,computing_panda = computing, precision=precision, alpha = alpha, th_motifs=th_motifs, delta=delta, tune_delta=tune_delta)
    print('All done!')

