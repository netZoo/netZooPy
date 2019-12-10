import pytest
import os
from netZooPy import condor
import pandas as pd
import numpy as np

def test_condor():
    print('Start Condor run ...')
    network        ='tutorials/condor/toynetwork.csv'
    output_file    ='tar_memb.txt'
    gt_file        ='tests/condor/travis_tar_memb.txt'

    #1. Vanilla panda
    condor.condor(network)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
