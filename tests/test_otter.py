import pytest
import os
from netZooPy import otter
from netZooPy.lioness.lioness_for_otter import *
from netZooPy.panda import Panda
import pandas as pd

import logging

LOGGER = logging.getLogger(__name__)

def test_otter():
    
    LOGGER.warning('Test Otter calculation')
    print("Start Otter run ...")
    W = "tests/otter/w.csv"
    W = pd.read_csv(W, header=None)
    W = W.values
    C = "tests/otter/c.csv"
    C = pd.read_csv(C, header=None)
    C = C.values
    P = "tests/otter/p.csv"
    P = pd.read_csv(P, header=None)
    P = P.values
    gt_file = "tests/otter/test_otter.csv"

    # 1. Call Otter
    W = otter.otter(W, P, C, Iter=1, lam=0.0035, gamma=0.335)
    gt = pd.read_csv(gt_file, header=None)
    W = pd.DataFrame(data=W)
    pd.testing.assert_frame_equal(W, gt, rtol=1e-10, check_exact=False)
    
    
    LOGGER.warning('Test Otter class')
    
    
    expression_fn = 'tests/puma/ToyData/ToyExpressionData.txt'
    ppi_fn = 'tests/puma/ToyData/ToyPPIData.txt'
    motif_fn = 'tests/puma/ToyData/ToyMotifData.txt'
    output_otter = './otter_test.txt'
    output_otter = './otter_test_panda.txt'

    lioobj = LionessOtter(expression_fn, motif_fn, ppi_fn, mode_process='intersection')
    lioobj.run_otter(output_otter)
    otter_otter = lioobj.all_otter
    
    # 1. Intersection
    panda_obj = Panda(
        expression_fn,
        motif_fn,
        ppi_fn,
        save_tmp=False,
        remove_missing=False,
        keep_expression_matrix=True,
        modeProcess="intersection",
    )
   
    W = panda_obj.motif_matrix
    P = panda_obj.ppi_matrix
    C = panda_obj.correlation_matrix
    
    panda_otter = otter.otter(W, P, C, Iter=60, lam=0.0035, gamma=0.335, eta = 0.00001, bexp = 1)
    
    
