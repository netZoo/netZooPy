import pytest
import os
from netZooPy.panda.panda import Panda
import pandas as pd

def test_panda():
    #print(os.getcwd())
    print('Start Panda run ...')
    ppi            ='tests/ToyData/ToyPPIData.txt'
    motif          ='tests/ToyData/ToyMotifData.txt'
    expression_data='tests/ToyData/ToyExpressionData.txt'
    lioness_file   =''
    rm_missing     = False
    output_file    ='travis_test_panda.txt'
    gt_file        ='tests/panda/test_panda.txt'
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file))
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(gt_file, sep='\t', header=None)
    gt =pd.read_csv(output_file, sep='\t', header=None)
    assert(gt.equals(round(res,6)))
    #pd.testing.assert_frame_equal(res,gt,check_less_precise=2,check_exact=False)
    print('Test panda passed was successful!')
