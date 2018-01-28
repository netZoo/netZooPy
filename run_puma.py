#!/usr/bin/env python
"""Run PANDA algorithm from the command line.
"""
import sys
import getopt
import pypanda

def main(argv):
    """ Run pypanda
    -h help
    -e (required) expression values
    -m (required) pair file of motif edges
    -p (required) pair file of PPI edges
    -o (required) output file
    Example:
    python run_puma -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o test_puma.txt
    """
    #create variables
    expression_data = None
    motif = None
    ppi = None
    mir = None
    output_file = None
    help_text = 'pypanda options:\n\
                \t-e, --expression (required) <expression_data.txt>\n\
                \t-m, --motif (required) <motif_data.txt>\n\
                \t-p (required) <ppi_data.txt>\n\
                \t-i, --mir (optional)<MiRList>\n\
                \t-o, --output (required) <output_file_name>'

    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:i:o:', ['help', 'expression=', 'motif=', 'ppi=', 'mir=', 'out='])
    except getopt.GetoptError:
        print(help_text)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(help_text)
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
    #check if required options are given
    if expression_data and motif and ppi and miR:
        print('Input data:')
        print('Expression:', expression_data)
        print('Motif data:', motif)
        print('PPI data:', ppi)
        print('miR data:', miR)
    else:
        print('Missing input file!')
        print(help_text)
        sys.exit()

    # Run puma
    print('Start Puma run ...')
    p = pypanda.Puma(expression_data, motif, ppi, miR, save_tmp=True)
    p.save_panda_results(output_file)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
