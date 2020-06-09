#!/usr/bin/env python
"""
Description:
  Run predScore after benchmark.

Usage:
  -h, --help: help
  -i, --indir: directory where intersections are location
  -o, --outdir
  -c, --cell: if specific cell line wanted
  -t, --TF: if specific TF wanted
  -o, --out: output folder
  -f, --format: output format (txt, npy, or mat)
  start: to start from nth sample (optional)
  end: to end at nth sample (optional, must with start)
  
Example:
  source /proj/relibs/relib00/conda/bin/activate
  source activate mypy3
  python run_predScore.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -o /tmp -f npy 1 2
"""


import sys
import getopt

def main(argv):
    #Create variables
    indir = None
    outdir = None
    cell = None
    TF = None
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:n:o:f:', ['help', 'indir=','outdir=','cell=','TF='])
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        print(__doc__)
        return 2

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            return 0
        elif opt in ('-i', '--indir'):
            indir = arg
        elif opt in ('-o', '--outdir'):
            outdir = arg
        elif opt in ('-c', '--cell'):
            cell = arg
        elif opt in ('-t','--TF'):
            TF = arg

        else:
            print('Unknown option', opt)
            return 1

    #Check if required options are given
    if indir is None or outdir is None:
        print('Missing argument!')
        print(__doc__)
        return 1
    else:
        print('indir: ' indir)
        print('outdir: ', outdir)
    if TF is not None and cell is not None:
        print('TF:   ', TF)
        print('cell:   ', cell)
    elif TF is not None and cell is None:
        print('TF:   ', TF)
    elif TF is None and cell is not None:
        print('cell:   ', cell)
    else:
        print('all cell and TF combos')

    # Run panda
    print('Start predScore ...')
    predScore(indir,outdir,cell,TF)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
