#!/usr/bin/env python3
"""
Description:
   Run MILIPEED algorithm from the command line.

Usage:
  -h,  --help: help
  -e,  --expression: expression values
  -me, --methylation: beta values
  -p,  --ppi: pair file of PPI edges
  -cg, --CGmap: listing TF - gene - CG in that order
  -o,  --out: output file
  start: to start from nth sample (optional)
  end: to end at nth sample (optional, must with start)
  
Example:
  python run_milipeed.py -e ./ToyData/ToyExpressionData.txt -me ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o test_milipeed.txt
"""
import sys
import getopt
from netZooPy.milipeed.milipeed import Milipeed
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
        opts, args = getopt.getopt(argv, 'he:m:p:o:rq:', ['help', 'expression=', 'methylation=', 'ppi=','cgmap=','out='])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            sys.exit()
        elif opt in ('-e', '--expression'):
            expression_data = arg
        elif opt in ('-me', '--methylation'):
            methyl = arg
        elif opt in ('-cg', '--cgmap'):
            cgmap = arg
        elif opt in ('-p', '--ppi'):
            ppi = arg
        elif opt in ('-o', '--out'):
            output_file = arg

    #Check if required options are given
    print('Input data:')
    print('Expression:', expression_data)
    print('Methylation:', methylation_data)
    print('PPI data:', ppi)
    print('CGMap data:', cgmap)
    print('Sample range: ', start, '-', end)

    if methylation_data is None or motif is None or expression_data is None:
        print('Missing methylation inputs!')
        print(__doc__)
        sys.exit()

    # Run MILIPEED
    print('Start Milipeed run ...')
    panda_obj = Milipeed(expression_data,methylation, motif, ppi, cgmap, outdir,out,start=start, end=end)

    panda_obj.save_milipeed_results(output_file)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



