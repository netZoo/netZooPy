import pytest
import os
from netZooPy.puma.puma import Puma
import pandas as pd
import numpy as np

def test_puma():
    #print(os.getcwd())
    print('Start Puma run ...')
    ppi            ='tests/puma/ToyData/ToyPPIData.txt'
    motif          ='tests/puma/ToyData/ToyMotifData.txt'
    expression_data='tests/puma/ToyData/ToyExpressionData.txt'
    mir_file       ='tests/puma/ToyData/ToyMiRList.txt'
    lioness_file   =''
    rm_missing     = False
    output_file    ='travis_test_panda.txt'
    gt_file        ='tests/puma/test_puma.txt'

    #1. Vanilla panda
    panda_obj      = Puma(expression_data, motif, ppi, mir_file,save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='legacy')
    panda_obj.save_puma_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #2. with argument values
    rm_missing= False
    print('Test puma was successful!')
