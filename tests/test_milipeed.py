import pytest
import os
from netZooPy.milipeed.milipeed import Milipeed
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np

def test_milipeed():
    print('Start milipeed test!')
    #1. First generate temporary PANDA files as inputs for Milipeed
    ppi_file             ='tests/puma/ToyData/ToyPPIData.txt'
    motif_file           ='tests/puma/ToyData/ToyMotifData.txt'
    methylation_file     ='test/milipeed/ToyData/ToyMethylationData.txt'
    expression_file      ='tests/milipeed/ToyData/ToyExpressionData.txt'
    milipeed_file        ='Travis'
    milipeed_Nocorr_file = 'Travis'
    rm_missing           =False
    output_file          ='milipeed.npy'
    milipeed_obj         = Panda(expression_file, motif_file, methylation_file, ppi_file, start=1,end=1)
    # Set parameters
    milipeed_obj = Milipeed(milipeed_obj)
    milipeed_obj.export_milipeed_results(milipeed_obj,milipeed_file)
    # Read first milipeed network
    res  = np.load('milipeed_output/milipeed.1.npy')
    gt = np.load('tests/lioness/lioness.1.npy')
    # Compare to ground truth
    assert(np.allclose(gt,res))

    #2. Testing Milipeed with exoression set to None to compute Milipeed on coexpression networks
    expression_data              = None
    milipeed_obj.expression_data = None
    # Make sure to keep epxression matrix for next step
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, modeProcess='legacy')
    lioness_obj    = Lioness(panda_obj)
    lioness_obj.save_lioness_results(lioness_file)

    milipeed_obj2      = Milipeed(milipeed_obj)
    milipeed_obj2.save_milipeed_results(milipeed_obj2,milipeed_Nocorr_file)
    # Read first milipeed network
    res  = np.load('milipeed_output/milipeed_noCorr.1.npy')

    gt   = np.load('tests/lioness/lionessCoexpression.1.npy')
    # Compare to ground truth
    assert(np.allclose(gt,res))
