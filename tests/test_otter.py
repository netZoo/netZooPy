import pytest
import os
from netZooPy import otter
import pandas as pd
import numpy as np

def test_otter():
    print('Start Otter run ...')
    W       = 'tests/otter/w.csv'
    W       = pd.read_csv(W, header=None)
    W       = W.values
    C       = 'tests/otter/c.csv'
    C       = pd.read_csv(C, header=None)
    C       = C.values
    P       = 'tests/otter/p.csv'
    P       = pd.read_csv(P, header=None)
    P       = P.values
    gt_file = 'tests/otter/test_otter.csv'

    #1. Call Otter
    W  =otter.otter(W,P,C,Iter=1, lam = 0.0035, gamma = 0.335)
    gt =pd.read_csv(gt_file, header=None)
    W  =pd.DataFrame(data=W)
    pd.testing.assert_frame_equal(W,gt,check_less_precise=10,check_exact=False)
