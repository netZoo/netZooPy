import pytest
import os
from netZooPy.lioness.lioness import Lioness
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np

def test_lioness():
    print('Start lioness test!')
    #1. First generate temporary PANDA files as inputs for Lioness
    ppi            ='tests/ToyData/ToyPPIData.txt'
    motif          ='tests/ToyData/ToyMotifData.txt'
    expression_data='tests/ToyData/ToyExpressionData.txt'
    lioness_file   ='Travis'
    rm_missing     =False
    output_file    ='panda.npy'
    gt_file        ='tests/panda/test_panda.txt'
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file))
    

    # Set parameters
    lioness_obj = Lioness(panda_obj)
    lioness_obj.save_lioness_results(lioness_file)

    # Read first lioness network
    gt  = np.load('lioness_output/lioness.1.npy')
    res = np.load('tests/lioness/lioness.1.npy')

    # Compare to ground truth
    assert(np.array_equal(gt,res))
