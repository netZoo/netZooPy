#!/usr/bin/env python
"""Run PUMA algorithm from the command line.

Usage:
run_puma
  -h, --help: help
  -e, --expression: expression values
  -m, --motif: pair file of motif edges
  -p, --ppi: pair file of PPI edges
  -i, --mir (required): miR file
  -o, --out: output file
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
    rm_missing = False
    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:i:o:r', ['help', 'expression=', 'motif=', 'ppi=', 'mir=', 'out=', 'rm_missing'])
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
    #Check if required options are given
    print('Input data:')
    print('Expression:', expression_data)
    print('Motif data:', motif)
    print('PPI data:', ppi)
    if not expression_data and not motif and not ppi and not miR:
        print('Missing inputs!')
        print(__doc__)
        sys.exit()

    # Run PUMA
    print('Start Puma run ...')
    puma_obj = pypanda.Puma(expression_data, motif, ppi, miR, save_tmp=True, remove_missing=rm_missing)
    puma_obj.save_puma_results(output_file)
    puma_obj.top_network_plot(top=100, file='puma_top100genes.png')
    #indegree = puma_obj.return_panda_indegree()
    #outdegree = puma_obj.return_panda_outdegree()
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



