#!/usr/bin/env python
"""Run PUMA algorithm from the command line.

Usage:
run_puma
  -h, --help: help
  -e, --expression (required): expression values
  -m, --motif (required): pair file of motif edges
  -p, --ppi (required): pair file of PPI edges
  -i, --mir (required): miR file
  -o, --out (required): output file
  Example:
  python run_puma.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o test_puma.txt
"""
import sys
import getopt
import pypanda

def main(argv):
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    miR = None
    output_file = "output_puma.txt"
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:i:o:', ['help', 'expression=', 'motif=', 'ppi=', 'mir=', 'out='])
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
    #Check if required options are given
    if expression_data and motif and ppi and miR:
        print('Input data:')
        print('Expression:', expression_data)
        print('Motif data:', motif)
        print('PPI data:', ppi)
        print('miR data:', miR)
    else:
        print('Missing inputs!')
        print(__doc__)
        sys.exit()

    # Run puma
    print('Start Puma run ...')
    p = pypanda.Puma(expression_data, motif, ppi, miR, save_tmp=True)
    p.save_puma_results(output_file)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
