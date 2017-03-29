#!/usr/bin/env python
"""Run PANDA algorithm from the command line.
"""
import sys
import getopt
import pypanda

def main(argv):
    """ Run pypanda.
    -h help
    -e (required) expression values
    -m (required) pair file of motif edges
    -p (required) pair file of PPI edges
    -o (required) output file
    """
    #create variables
    expression_data = None
    motif = None
    ppi = None
    output_file = None
    help_text = 'pypanda options:\n\
                \t-e (required) <expression_data.txt>\n\
                \t-m (required) <motif_data.txt>\n\
                \t-p (required) <ppi_data.txt>\n\
                \t-o (required) <output_file_name>'

    # Get input options
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:o:', ['help', 'expression=', 'motif=', 'ppi=', 'out='])
    except getopt.GetoptError:
        print help_text
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print help_text
            sys.exit()
        if opt == '-e':
            expression_data = arg
        if opt == '-m':
            motif = arg
        if opt == '-p':
            ppi = arg
        if opt == '-o':
            output_file = arg

    #check if required options are given
    if expression_data and motif and ppi:
        print 'Input data:'
        print 'Expression:', expression_data
        print 'Motif data:', motif
        print 'PPI data:', ppi
    else:
        print 'Missing input file!'
        print help_text
        sys.exit()

    # Run panda
    print 'Start Panda run ...'
    p = pypanda.Panda(expression_data, motif, ppi, save_tmp=True)
    p.save_panda_results(output_file)
    print 'All done!'

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
