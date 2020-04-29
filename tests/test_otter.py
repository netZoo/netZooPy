import pytest
import os
from netZooPy import otter
import pandas as pd
import numpy as np

def test_condor():
    print('Start Otter run ...')
    W              = 'tests/otter/W.csv'
    C              = 'tests/otter/C.csv'
    P              = 'tests/otter/P.csv'
    gt_file        = 'tests/otter/test_otter.csv'

    #1. Vanilla panda
    W=otter.otter(W,C,P)
    gt =pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(W,gt,check_less_precise=False,check_exact=False)
