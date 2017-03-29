#!/usr/bin/env python
"""Run LIONESS algorithm from the command line.

Usage:

./run_lioness.py -e <expression_matrix> -m <normalized_motif_matrix> -p <normalized_ppi_matrix> -n <panda_network> -o <output_dir> -f <output_fmt> [start end]

  -e: expression matrix (.npy)
  -m: motif matrix, normalized (.npy)
  -p: ppi matrix, normalized (.npy)
  -n: PANDA network (.npy)
  -o: output folder
  -f: output format (txt, npy, or mat)
  start: to start from nth sample (optional)
  end: to end at nth sample (optional, must with start)

E.g., ./run_lioness.py -e expression.npy -m motif.npy -p ppi.npy -n panda.npy -o /tmp -f npy 1 100
"""
import sys
import getopt
import pypanda

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "he:m:p:n:o:f:")
    except getopt.GetoptError as err:
        print str(err)  # will print something like "option -a not recognized"
        print __doc__
        return 2

    for x,y in opts:
        if x == '-h':
            print __doc__
            return 0
        elif x == '-e':
            expression = y
        elif x == '-m':
            motif = y
        elif x == '-p':
            ppi = y
        elif x == '-n':
            panda_net = y
        elif x == '-o':
            save_dir = y
        elif x == '-f':
            save_fmt = y
        else:
            print 'Unknown option', x
            return 1

    if len(args) == 2:
        start, end = map(int, args)
    else:
        start, end = 1, None

    try:
        print 'Input data:'
        print 'Expression:   ', expression
        print 'Motif matrix: ', motif
        print 'PPI matrix:   ', ppi
        print 'PANDA network:', panda_net
        print 'Output folder:', save_dir
        print 'Output format:', save_fmt
        print 'Sample range: ', start, '-', end
    except UnboundLocalError as err:
        print 'Missing argument!'
        print __doc__
        return 1

    # Run panda
    print 'Start LIONESS run ...'
    L = pypanda.Lioness(expression, motif, ppi, panda_net, start=start, end=end, save_dir=save_dir, save_fmt=save_fmt)
    print 'All done!'

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
