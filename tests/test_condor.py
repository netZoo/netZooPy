import pytest
import os
from netZooPy import condor
import pandas as pd
import numpy as np

def test_condor():
    print('Start Condor run ...')
    network        ='tutorials/condor/toynetwork.csv'
    output_file1   ='tar_memb.txt'
    output_file2   ='reg_memb.txt'
    gt_file1       ='tests/condor/gh_tar_memb.txt'
    gt_file2       ='tests/condot/gh_reg_memb.txt'

    #1. Condor
    condor.run_condor(network)
    res1=pd.read_csv(output_file1, sep=' ', header=None)
    res2=pd.read_csv(output_file2, sep=' ', header=None)
    gt1 =pd.read_csv(gt_file1, sep=' ', header=None)
    gt2 =pd.read_csv(gt_file2, sep=' ', header=None)
    pd.testing.assert_frame_equal(res1,gt1,check_less_precise=False,check_exact=False)
    pd.testing.assert_frame_equal(res2,gt2,check_less_precise=False,check_exact=False)
