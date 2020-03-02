#!/usr/bin/env python
"""
Description:
   Run MILIPEED algorithm from the command line.

Usage:
  -h, --help: help
  -e, --expression: expression values
  -m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided
  -me,--methylation: beta values
  -p, --ppi: pair file of PPI edges
  -o, --out: output file
  -r, --rm_missing
  start: to start from nth sample (optional)
  end: to end at nth sample (optional, must with start)
  
Example:
  python run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o test_panda.txt -q output_panda.txt
"""
import sys
import getopt
from netZooPy.netZooPy.milipeed import Milipeed
from netZooPy.panda.panda import Panda


def main(argv):
    #Create variables
    expression_data = None
    motif = None
    methylation_data = None
    ppi = None
    output_file = "output_milipeed.txt"
    rm_missing = False
    lioness_file = False
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:o:rq:', ['help', 'expression=', 'motif=', 'methylation=', 'ppi=', 'out='])
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
        elif opt in ('-me', '--methylation'):
            methyl = arg
        elif opt in ('-p', '--ppi'):
            ppi = arg
        elif opt in ('-o', '--out'):
            output_file = arg
        elif opt in ('-n'):
            panda_net = arg
    #Check if required options are given
    print('Input data:')
    print('Expression:', expression_data)
    print('Motif data:', motif)
    print('Methylation:', methylation_data)
    print('PPI data:', ppi)
    print('Sample range: ', start, '-', end)

    if methylation_data is None or motif is None or expression_data is None:
        print('Missing methylation inputs!')
        print(__doc__)
        sys.exit()

    # Run PANDA
    print('Start Milipeed run ...')
    panda_obj = Milipeed(expression_data,methylation, motif, ppi, save_tmp=True, start=start, end=end)
    #panda_obj = pypanda.Panda(expression_data, motif, None, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(None, motif, ppi, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(None, motif, None, save_tmp=True, remove_missing=rm_missing)
    #panda_obj = pypanda.Panda(expression_data, None, ppi, save_tmp=True, remove_missing=rm_missing)
    panda_obj.save_milipeed_results(output_file)
    #panda_obj.top_network_plot(top=70, file='panda_topgenes.png')
    #indegree = panda_obj.return_panda_indegree()
    #outdegree = panda_obj.return_panda_outdegree()

    # if lioness_file:
    #     from netZooPy.lioness.lioness import Lioness
    #     lioness_obj = Lioness(panda_obj)
    #     lioness_obj.save_lioness_results(lioness_file)
    # print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



