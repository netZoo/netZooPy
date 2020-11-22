#!/usr/bin/env python

import sys
import getopt
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda

def main(argv):
    """
    Description:
        Run LIONESS algorithm from the command line.

    Usage:
        -h, --help: help
        -e, --expression: expression matrix (.npy)
        -m, --motif: motif matrix, normalized (.npy)
        -p, --ppi: ppi matrix, normalized (.npy)
        -g, --comp: use cpu (default) or gpu
        -r, --pre: number of digits to calcluate
        -c, --ncores: number cores
        -n, --npy: PANDA network (.npy)
        -o, --out: output folder
        -f, --format: output format (txt, npy, or mat)
        start: to start from nth sample (optional)
        end: to end at nth sample (optional, must with start)
    
    Example:
        python3 run_lioness.py -e ../../tests/ToyData/ToyExpressionData.txt -m ../../tests/ToyData/ToyMotifData.txt -p ../../tests/ToyData/ToyPPIData.txt -g cpu -r single -c 2 -o /tmp -f npy 1 2

    Reference:
        Kuijjer, Marieke Lydia, et al. "Estimating sample-specific regulatory networks." Iscience 14 (2019): 226-240.
    """
    #Create variables
    expression_data = None
    motif = None
    ppi = None
    comp = None
    pre = None
    ncores = None
    save_dir = None
    save_fmt = None
    try:
        opts, args = getopt.getopt(argv, 'he:m:p:g:r:c:n:o:f:', ['help', 'expression=', 'motif=','ppi=','comp=','pre=','ncores=', 'out=', 'format='])
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
        elif opt in ('-g', '--comp'):
            comp = arg
        elif opt in ('-r', '--pre'):
            pre = arg
        elif opt in ('-c', '--ncores'):
            ncores = arg
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
        print('compute core: ', comp)
        print('precision:    ', pre)
        print('n cores:      ', ncores)        
        print('Output folder:', save_dir)
        print('Output format:', save_fmt)
        print('Sample range: ', start, '-', end)

    # Run panda
    print('Start LIONESS run ...')
    obj = Panda(expression_data, motif, ppi, keep_expression_matrix=True,save_memory=False)
    L   = Lioness(obj, computing=comp, precision=pre,ncores=ncores,start=start, end=end, save_dir=save_dir, save_fmt=save_fmt)
    print('All done!')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
