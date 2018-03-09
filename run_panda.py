#!/usr/bin/env python
"""Run PANDA algorithm from the command line.

Usage:
  -h, --help: help
  -e, --expression (required): expression values
  -m, --motif (required): pair file of motif edges
  -p, --ppi (required): pair file of PPI edges
  -o, --out (required): output file
"""
import sys
import getopt
import pypanda

def main(argv):
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    output_file = "output_panda.txt"
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:o:', ['help', 'expression=', 'motif=', 'ppi=', 'out='])
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

    #Check if required options are given
    if expression_data and motif and ppi:
        print('Input data:')
        print('Expression:', expression_data)
        print('Motif data:', motif)
        print('PPI data:', ppi)
    else:
        print('Missing input file!')
        print(__doc__)
        sys.exit()

    # Run panda
    print('Start Panda run ...')
    p = pypanda.Panda(expression_data, motif, ppi, save_tmp=True)
    p.save_panda_results(output_file)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
