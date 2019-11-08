#!/usr/bin/env python
"""
Description:
  Run LIONESS algorithm from the command line.

Usage:
  -h, --help: help
  -e, --expression: expression matrix (.npy)
  -m, --motif: motif matrix, normalized (.npy)
  -p, --ppi: ppi matrix, normalized (.npy)
  -n, --npy: PANDA network (.npy)
  -o, --out: output folder
  -f, --format: output format (txt, npy, or mat)
  start: to start from nth sample (optional)
  end: to end at nth sample (optional, must with start)
  
Example:
  python3 run_lioness.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -o /tmp -f npy 1 2
"""
import sys
import getopt
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda

def main(argv):
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    save_dir = None
    save_fmt = None
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:n:o:f:', ['help', 'expression=', 'motif=', 'ppi=', 'out=', 'format='])
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        print(__doc__)
        return 2

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            return 0
        elif opt in ('-e', '--expression'):
            expression_data = arg
        elif opt in ('-m', '--motif'):
            motif = arg
        elif opt in ('-p', '--ppi'):
            ppi = arg
        elif opt in ('-n'):
            panda_net = arg
        elif opt in ('-o', '--out'):
            save_dir = arg
        elif opt in ('-f', '--format'):
            save_fmt = arg
        else:
            print('Unknown option', opt)
            return 1

    start, end = 1, None
    if len(args) == 2:
        start, end = map(int, args)

    #Check if required options are given
    if expression_data is None or motif is None or ppi is None \
            or save_dir is None or save_fmt is None:
        print('Missing argument!')
        print(__doc__)
        return 1
    else:
        print('Input data:')
        print('Expression:   ', expression_data)
        print('Motif matrix: ', motif)
        print('PPI matrix:   ', ppi)
        print('Output folder:', save_dir)
        print('Output format:', save_fmt)
        print('Sample range: ', start, '-', end)

    # Run panda
    print('Start LIONESS run ...')
    obj = Panda(expression_data, motif, ppi, keep_expression_matrix=True)
    L   = Lioness(obj, start=start, end=end, save_dir=save_dir, save_fmt=save_fmt)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
