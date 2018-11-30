#!/usr/bin/env python
"""Run PUMA algorithm from the command line.

Usage:
run_puma
  -h, --help: help
  -e, --expression: expression values
  -m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided
  -p, --ppi: pair file of PPI edges
  -i, --mir (required): miR file
  -o, --out: output file
  -r, --rm_missing
  -q, --lioness: output for Lioness single sample networks 
  Example:
  python run_puma.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o test_puma.txt -q output_lioness.txt
"""
import sys
import getopt
import pypuma

def main(argv):
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    miR = None
    output_file = "output_puma.txt"
    rm_missing = False
    lioness_file = False
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:i:o:rq:', ['help', 'expression=', 'motif=', 'ppi=', 'mir=', 'out=', 'rm_missing', 'lioness'])
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
        elif opt in ('-i', '--mir'):
            miR = arg
        elif opt in ('-r', '--rm_missing'):
            rm_missing = arg
        elif opt in ('-q', '--lioness'):
            lioness_file = arg
    #Check if required options are given
    print('Input data:')
    print('Expression:', expression_data)
    print('Motif data:', motif)
    print('PPI data:', ppi)
    if (expression_data is None and motif is None) or miR is None:
        print('Missing inputs!')
        print(__doc__)
        sys.exit()

    # Run PUMA
    print('Start Puma run ...')
    puma_obj = pypuma.Puma(expression_data, motif, ppi, miR, save_tmp=True, remove_missing=rm_missing, keep_expression_matrix=bool(lioness_file))
    puma_obj.save_puma_results(output_file)
    #puma_obj.top_network_plot(top=100, file='puma_top100genes.png')
    #indegree = puma_obj.return_panda_indegree()
    #outdegree = puma_obj.return_panda_outdegree()

    if lioness_file:
        from pypuma.lioness import Lioness
        lioness_obj = Lioness(puma_obj)
        lioness_obj.save_lioness_results(lioness_file)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



