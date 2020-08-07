#!/usr/bin/env python

import sys
import getopt
from netZooPy.panda.panda import Panda

def main(argv):
    """
    Description:
        Run PANDA algorithm from the command line.

    Inputs:
        -h, --help: help
        -e, --expression: expression values
        -m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided
        -p, --ppi: pair file of PPI edges
        -o, --out: output file
        -r, --rm_missing
        -q, --lioness: output for Lioness single sample networks 
    
    Example:
        python run_panda.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -o test_panda.txt -q output_panda.txt

    Reference:
        Glass, Kimberly, et al. "Passing messages between biological networks to refine predicted interactions." PloS one 8.5 (2013): e64832.
    """
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    output_file = "output_panda.txt"
    rm_missing = False
    lioness_file = False
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:o:rq:', ['help', 'expression=', 'motif=', 'ppi=', 'out=', 'rm_missing', 'lioness'])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            sys.exit()
        elif opt in ('-e', '--expression'):
            expression_data = arg
        elif opt in ('-m', '--motif'):
            motif = arg
        elif opt in ('-p', '--ppi'):
            ppi = arg
        elif opt in ('-o', '--out'):
            output_file = arg
        elif opt in ('-r', '--rm_missing'):
            rm_missing = arg
        elif opt in ('-q', '--lioness'):
            lioness_file = arg
    #Check if required options are given
    print('Input data:')
    print('Expression:', expression_data)
    print('Motif data:', motif)
    print('PPI data:', ppi)
    if expression_data is None and motif is None:
        print('Missing inputs!')
        print(__doc__)
        sys.exit()

    # Run PANDA
    print('Start Panda run ...')
    panda_obj = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing, keep_expression_matrix=bool(lioness_file))
    #panda_obj = pypanda.Panda(expression_data, motif, None, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(None, motif, ppi, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(None, motif, None, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(expression_data, None, ppi, save_tmp=True, remove_missing=rm_missing)
    panda_obj.save_panda_results(output_file)
    #panda_obj.top_network_plot(top=70, file='panda_topgenes.png')
    #indegree = panda_obj.return_panda_indegree()
    #outdegree = panda_obj.return_panda_outdegree()

    if lioness_file:
        from netZooPy.lioness.lioness import Lioness
        lioness_obj = Lioness(panda_obj)
        lioness_obj.save_lioness_results(lioness_file)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))